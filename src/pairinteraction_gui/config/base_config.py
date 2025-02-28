from typing import TYPE_CHECKING

from pairinteraction_gui.qobjects import WidgetV

if TYPE_CHECKING:
    from pairinteraction_gui.page.base_page import SimulationPage


class BaseConfig(WidgetV):
    """Base class for configurations."""

    title: str

    def __init__(self, parent: "SimulationPage") -> None:
        """Initialize the base section."""
        self.page = parent
        super().__init__(parent)
