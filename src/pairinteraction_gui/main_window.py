import logging
from typing import TYPE_CHECKING, Optional, TypeVar

from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QAction, QActionGroup, QIcon, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDockWidget,
    QMainWindow,
    QSizePolicy,
    QStatusBar,
    QToolBar,
    QWidget,
)

from pairinteraction_gui.page import (
    AboutPage,
    LifetimesPage,
    SystemAtomPage,
    SystemPairPage,
)
from pairinteraction_gui.qobjects import NamedStackedWidget

if TYPE_CHECKING:
    from pairinteraction_gui.page.base_page import BasePage

    ChildType = TypeVar("ChildType")

logger = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    """Main window for the PairInteraction GUI application."""

    def __init__(self) -> None:
        """Initialize the main window."""
        super().__init__()

        self.setWindowTitle("PairInteraction")
        self.resize(1200, 800)

        self.apply_modern_style()

        self.statusbar = self.setup_statusbar()
        self.dockwidget = self.setup_dockwidget()

        self.stacked_pages = self.setup_stacked_pages()
        self.toolbar = self.setup_toolbar()

        self.init_keyboard_shortcuts()

    def findChild(  # type: ignore
        self, type_: type["ChildType"], name: str, options: Optional["Qt.FindChildOption"] = None
    ) -> "ChildType":  # overwrite type hints
        if options is None:
            options = Qt.FindChildOption.FindChildrenRecursively
        return super().findChild(type_, name, options)  # type: ignore

    def apply_modern_style(self) -> None:
        """Apply modern styling to the application."""
        self.setStyleSheet("""
            QMainWindow {
                background-color: #f8f9fa;
            }
            QLabel {
                color: #212529;
                font-family: 'Segoe UI', 'Arial', sans-serif;
            }
            QStatusBar {
                background-color: #343a40;
                color: #f8f9fa;
                font-size: 12px;
                padding: 3px;
            }
            QDockWidget {
                font-family: 'Segoe UI', 'Arial', sans-serif;
                font-weight: 500;
                color: #212529;
            }
            QDockWidget::title {
                background-color: #e9ecef;
                padding: 6px;
                border: none;
                border-bottom: 1px solid #dee2e6;
            }
            QToolBar {
                border: none;
                spacing: 15px;
            }
            QToolButton {
                background-color: transparent;
            }
        """)

    def setup_statusbar(self) -> QStatusBar:
        """Set up the status bar.

        The status bar message is set to "Ready" by default.
        It can be updated with a new message by either from the main window instance:
            `self.statusbar.showMessage("Ready", timeout=0)`
        or from outside the main window instance:
            `QApplication.sendEvent(self, QStatusTipEvent("Ready"))`
        """
        statusbar = QStatusBar(self)
        statusbar.setFixedHeight(25)
        self.setStatusBar(statusbar)
        statusbar.showMessage("Ready", timeout=0)
        return statusbar

    def setup_dockwidget(self) -> QDockWidget:
        """Create a configuration dock widget for the main window."""
        dockwidget = QDockWidget()
        dockwidget.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea)
        dockwidget.setTitleBarWidget(QWidget())  # This removes the title bar

        dockwidget.setMinimumWidth(350)
        dockwidget.setStyleSheet("""
            QToolBox {
                background-color: white;
                border: 1px solid #dee2e6;
                border-radius: 5px;
            }
            QToolBox::tab {
                background-color: #343a40;
                border: 1px solid #dee2e6;
                border-radius: 4px;
                color: #f8f9fa;
                font-weight: bold;
                font-size: 15px;
            }
            QToolBox::tab:selected {
                background-color: #007bff;
            }
            QToolBox::tab:hover:!selected {
                background-color: #495057;
            }
            QLabel {
                color: black;
                font-size: 14px;
            }
        """)
        dockwidget.setVisible(False)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dockwidget)
        return dockwidget

    def setup_stacked_pages(self) -> NamedStackedWidget["BasePage"]:
        """Set up the different pages for each toolbar option."""
        stacked_pages = NamedStackedWidget["BasePage"]()
        self.setCentralWidget(stacked_pages)

        stacked_pages.addNamedWidget(SystemAtomPage(), "system_atom")
        stacked_pages.addNamedWidget(SystemPairPage(), "system_pair")
        stacked_pages.addNamedWidget(LifetimesPage(), "lifetimes")
        # stacked_pages.addNamedWidget(C6Page(), "c6")

        # stacked_pages.addNamedWidget(SettingsPage(), "settings")
        stacked_pages.addNamedWidget(AboutPage(), "about")
        return stacked_pages

    def setup_toolbar(self) -> QToolBar:
        """Set up the toolbar with icon buttons."""
        toolbar = QToolBar("Sidebar")
        toolbar.setMovable(False)
        toolbar.setOrientation(Qt.Orientation.Vertical)
        toolbar.setIconSize(QSize(32, 32))
        toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        toolbar.setStyleSheet("""
            QToolBar {
                background-color: #343a40;
                border: none;
                spacing: 15px;
                padding: 10px 5px;
            }
            QToolButton {
                border: none;
                padding: 8px;
                margin: 5px;
                color: #f8f9fa;
                font-weight: bold;
                font-size: 15px;
            }
            QToolButton:hover {
                background-color: #495057;
            }
            QToolButton:pressed {
                background-color: #212529;
            }
            QToolButton:checked {
                background-color: #007bff;
            }
        """)

        toolbar_group = QActionGroup(self)
        toolbar_group.setExclusive(True)

        for name, page in self.stacked_pages.items():
            # add a spacer widget
            if name == "about":
                spacer_widget = QWidget()
                spacer_widget.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
                toolbar.addWidget(spacer_widget)

            action = QAction(self)
            action.setObjectName(name)
            action.setText(page.title)
            action.setToolTip(page.tooltip)
            action.setCheckable(True)
            if page.icon_path:
                action.setIcon(QIcon(str(page.icon_path)))

            toolbar.addAction(action)
            toolbar_group.addAction(action)

            action.triggered.connect(lambda checked, name=name: self.stacked_pages.setCurrentNamedWidget(name))

        default_page = "system_atom"
        self.findChild(QAction, default_page).setChecked(True)
        self.stacked_pages.setCurrentNamedWidget(default_page)

        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, toolbar)

        return toolbar

    def init_keyboard_shortcuts(self) -> None:
        """Initialize keyboard shortcuts."""
        # Add Ctrl+W shortcut to close the window
        close_shortcut = QShortcut(QKeySequence("Ctrl+W"), self)
        close_shortcut.activated.connect(self.close)
