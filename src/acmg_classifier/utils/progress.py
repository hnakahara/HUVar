from __future__ import annotations
from contextlib import contextmanager
from typing import Generator, Iterable, TypeVar

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


def make_progress() -> Progress:
    return Progress(
        SpinnerColumn(),
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        MofNCompleteColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
    )


@contextmanager
def progress_iter(
    iterable: Iterable[T],
    description: str,
    total: int,
) -> Generator[Iterable[T], None, None]:
    """Wrap an iterable with a rich progress bar."""
    with make_progress() as progress:
        task = progress.add_task(description, total=total)

        def _tracked() -> Iterable[T]:
            for item in iterable:
                yield item
                progress.advance(task)

        yield _tracked()
