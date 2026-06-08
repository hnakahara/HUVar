"""Console progress-bar helpers (rich).

Two public APIs cover every loop in the codebase:

  * ``track(iterable, description, total=...)`` — drop-in iterator wrapper
    for simple ``for`` loops.
  * ``progress_bar(description, total=...)`` — context manager that yields
    an ``advance(n=1)`` callable, for cases where iteration is not in a
    plain loop (e.g. ``ThreadPoolExecutor.as_completed``).

Both auto-detect tty: when stderr is not a TTY (CI, redirected output)
they become no-ops so log files don't fill up with carriage-returns.
``set_enabled(True/False)`` overrides the auto-detection — used by the
``--no-progress`` CLI flag.
"""
from __future__ import annotations

import sys
from contextlib import contextmanager
from typing import Callable, Iterable, Iterator, Optional, TypeVar

from rich.console import Console
from rich.progress import (
    BarColumn,
    MofNCompleteColumn,
    Progress,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

T = TypeVar("T")

# Cap the redraw rate. rich defaults to 10 Hz AND repaints continuously to
# animate the spinner; profiling showed table/segment rendering eating ~10s
# over a fast 100-variant loop. 4 Hz is smooth enough for humans and roughly
# halves the render cost. (The spinner column is also dropped below for the
# same reason — the bar itself already shows liveness.)
_REFRESH_HZ = 4

# Module-level override set by the CLI (`--no-progress`) before the first
# progress-aware call. `None` keeps the automatic tty detection in effect.
_FORCE: Optional[bool] = None

# Number of live displays (bars/spinners) currently active. rich allows only
# one live display at a time, so `status()` consults this to skip its spinner
# when a bar is already on screen — a console-instance-agnostic check.
_LIVE_DEPTH = 0


def set_enabled(enabled: Optional[bool]) -> None:
    """Force progress bars on (True) / off (False), or restore auto (None)."""
    global _FORCE
    _FORCE = enabled


def _is_enabled() -> bool:
    """Return True if a progress bar should actually be drawn."""
    if _FORCE is not None:
        return _FORCE
    # stderr is used so that the bar does not corrupt a redirected stdout
    # (most pipelines pipe stdout to a file).
    return sys.stderr.isatty()


def is_enabled() -> bool:
    """Public predicate: will progress bars actually render this run?

    Lets callers skip extra work that only exists to feed a bar (e.g. polling a
    subprocess's output file) when bars are disabled (non-tty / --no-progress)."""
    return _is_enabled()


def make_progress() -> Progress:
    """Build a Progress instance with our standard column layout.

    Kept public so callers that need multiple concurrent tasks
    (e.g. one bar per chromosome) can manage a single Progress themselves.
    """
    return Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        # stderr keeps the bar out of redirected stdout (TSV output goes there).
        console=Console(stderr=True),
        # Leave the completed bar visible so the user sees the final timing.
        transient=False,
        # Throttle redraws (see _REFRESH_HZ) to cut rich's rendering overhead.
        refresh_per_second=_REFRESH_HZ,
    )


def track(
    iterable: Iterable[T],
    description: str,
    total: Optional[int] = None,
) -> Iterator[T]:
    """Iterate `iterable` with a progress bar (no-op when disabled).

    Prefers an explicit `total` because many bioinformatics iterables
    (file generators, tabix fetches) don't expose `len()`. Falls back to
    `len()` only when `total` is omitted AND the iterable is sized.
    """
    if not _is_enabled():
        yield from iterable
        return

    if total is None:
        try:
            total = len(iterable)  # type: ignore[arg-type]
        except TypeError:
            total = None

    global _LIVE_DEPTH
    with make_progress() as progress:
        task = progress.add_task(description, total=total)
        _LIVE_DEPTH += 1
        try:
            for item in iterable:
                yield item
                progress.advance(task)
        finally:
            _LIVE_DEPTH -= 1


@contextmanager
def status(description: str) -> Iterator[None]:
    """Indeterminate spinner for a single slow blocking call (no count).

    Use for opaque one-shot work like opening a multi-GB file, where there is
    nothing to count. No-op when progress is disabled, and silently skips the
    spinner if another live display (e.g. a progress bar) is already active —
    rich allows only one live display at a time, so nesting would error."""
    global _LIVE_DEPTH
    if not _is_enabled() or _LIVE_DEPTH > 0:
        # Disabled, or a bar is already on screen (only one live display is
        # allowed) — run the work without a spinner.
        yield
        return
    with Console(stderr=True).status(description):
        _LIVE_DEPTH += 1
        try:
            yield
        finally:
            _LIVE_DEPTH -= 1


@contextmanager
def progress_bar(
    description: str,
    total: Optional[int],
) -> Iterator[Callable[..., None]]:
    """Context manager yielding an ``advance(n=1)`` callable.

    Use when progress cannot be measured by simple iteration — typically a
    thread pool consumed via ``as_completed`` where the order of completion
    is not the order of submission. Pass ``total=None`` for an
    indeterminate bar (spinner + running count) when the total work is not
    known up-front, e.g. streaming XML parses.

    Example::

        with progress_bar("Annotating", total=len(variants)) as advance:
            for future in as_completed(futures):
                ...
                advance()
    """
    if not _is_enabled():
        # No-op advance keeps call sites identical between interactive and
        # non-interactive runs.
        yield lambda n=1: None
        return

    global _LIVE_DEPTH
    with make_progress() as progress:
        task = progress.add_task(description, total=total)
        _LIVE_DEPTH += 1
        try:
            yield lambda n=1: progress.advance(task, n)
        finally:
            _LIVE_DEPTH -= 1
