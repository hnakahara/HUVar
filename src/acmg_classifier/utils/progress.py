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
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

T = TypeVar("T")

# Module-level override set by the CLI (`--no-progress`) before the first
# progress-aware call. `None` keeps the automatic tty detection in effect.
_FORCE: Optional[bool] = None


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


def make_progress() -> Progress:
    """Build a Progress instance with our standard column layout.

    Kept public so callers that need multiple concurrent tasks
    (e.g. one bar per chromosome) can manage a single Progress themselves.
    """
    return Progress(
        SpinnerColumn(),
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

    with make_progress() as progress:
        task = progress.add_task(description, total=total)
        for item in iterable:
            yield item
            progress.advance(task)


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

    with make_progress() as progress:
        task = progress.add_task(description, total=total)
        yield lambda n=1: progress.advance(task, n)
