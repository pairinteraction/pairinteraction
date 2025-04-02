# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

from typing import TYPE_CHECKING, Generic, Optional, TypeVar

from PySide6.QtCore import QObject, Qt
from PySide6.QtWidgets import QFormLayout, QHBoxLayout, QLayout, QVBoxLayout, QWidget

if TYPE_CHECKING:
    ChildType = TypeVar("ChildType", bound=QObject)
    from pairinteraction_gui.main_window import MainWindow

LayoutType = TypeVar("LayoutType", bound=QLayout)


class Widget(QWidget, Generic[LayoutType]):
    """Custom Widget class."""

    layout_type: type[LayoutType]
    margin: tuple[int, int, int, int] = (0, 0, 0, 0)
    spacing: int = 0

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        *,
        name: Optional[str] = None,
        margin: Optional[tuple[int, int, int, int]] = None,
        spacing: Optional[int] = None,
    ) -> None:
        """Initialize the base section."""
        super().__init__(parent)

        if name is not None:
            self.setObjectName(name)
        if margin is not None:
            self.margin = margin
        if spacing is not None:
            self.spacing = spacing

        layout = self.layout_type(self)  # can be accessed as self.layout()
        layout.setContentsMargins(*self.margin)
        layout.setSpacing(self.spacing)

        self.setupWidget()
        self.postSetupWidget()

    def layout(self) -> LayoutType:
        return super().layout()  # type: ignore [return-value] # explicitly override type hints

    def window(self) -> "MainWindow":
        return super().window()  # type: ignore [return-value] # explicitly override type hints

    def findChild(  # type: ignore [override] # explicitly override type hints
        self, type_: type["ChildType"], name: str, options: Optional["Qt.FindChildOption"] = None
    ) -> "ChildType":
        if options is None:
            options = Qt.FindChildOption.FindChildrenRecursively
        return super().findChild(type_, name, options)  # type: ignore [return-value] # explicitly override type hints

    def setupWidget(self) -> None:
        """Set up the UI components.

        This method should be overwritten by subclasses to set up the UI components.
        """

    def postSetupWidget(self) -> None:
        """Post-process the UI components.

        This method should be overwritten by subclasses to post-process the UI components.
        """


class WidgetV(Widget[QVBoxLayout]):
    """Custom Widget class with vertical (QVBoxLayout) layout."""

    layout_type = QVBoxLayout


class WidgetH(Widget[QHBoxLayout]):
    """Custom Widget class with horizontal (QHBoxLayout) layout."""

    layout_type = QHBoxLayout


class WidgetForm(Widget[QFormLayout]):
    """Custom Widget class with form (QFormLayout) layout."""

    layout_type = QFormLayout
    margin = (5, 15, 5, 5)
    spacing = 10
