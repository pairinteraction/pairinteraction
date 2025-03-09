import logging
from collections.abc import ItemsView
from typing import Generic, Optional, TypeVar

from PySide6.QtWidgets import QStackedWidget, QWidget

logger = logging.getLogger(__name__)


WidgetType = TypeVar("WidgetType", bound=QWidget)


class NamedStackedWidget(QStackedWidget, Generic[WidgetType]):
    def __init__(self, parent: Optional[WidgetType] = None) -> None:
        super().__init__(parent)
        self._widgets: dict[str, WidgetType] = {}

    def addNamedWidget(self, widget: WidgetType, name: str) -> None:
        widget.setObjectName(name)
        self.addWidget(widget)
        self._widgets[name] = widget

    def setCurrentNamedWidget(self, name: str) -> None:
        widget = self.getNamedWidget(name)
        current_widget = self.currentWidget()
        if widget == current_widget:
            return
        self.setCurrentWidget(widget)
        logger.debug("Switched NamedStackedWidget to %s", name)

    def getNamedWidget(self, name: str) -> WidgetType:
        return self._widgets[name]

    def items(self) -> ItemsView[str, WidgetType]:
        return self._widgets.items()

    # just overwriting type hints
    def currentWidget(self) -> WidgetType:  # overwrite type hints
        return super().currentWidget()  # type: ignore

    def addWidget(self, widget: WidgetType) -> int:  # type: ignore  # overwrite type hints
        return super().addWidget(widget)
