# SPDX-FileCopyrightText: 2025 Pairinteraction Developers
# SPDX-License-Identifier: LGPL-3.0-or-later

import logging
from typing import TYPE_CHECKING, Optional, TypeVar

from PySide6.QtCore import QObject, QSize, Qt
from PySide6.QtGui import QAction, QActionGroup, QCloseEvent, QIcon, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QDockWidget,
    QMainWindow,
    QMessageBox,
    QSizePolicy,
    QStatusBar,
    QToolBar,
    QWidget,
)

from pairinteraction._wrapped import Database
from pairinteraction_gui.app import Application
from pairinteraction_gui.page import (
    LifetimesPage,
    OneAtomPage,
    TwoAtomsPage,
)
from pairinteraction_gui.page.base_page import SimulationPage
from pairinteraction_gui.qobjects import NamedStackedWidget
from pairinteraction_gui.theme import main_theme
from pairinteraction_gui.utils import download_databases_mp
from pairinteraction_gui.worker import MultiProcessWorker, MultiThreadWorker

if TYPE_CHECKING:
    from pairinteraction_gui.page import BasePage

    ChildType = TypeVar("ChildType", bound=QObject)

logger = logging.getLogger(__name__)


class MainWindow(QMainWindow):
    """Main window for the PairInteraction GUI application."""

    def __init__(self) -> None:
        """Initialize the main window."""
        super().__init__()

        self.setWindowTitle("PairInteraction")
        self.resize(1200, 800)
        self.setStyleSheet(main_theme)

        self.statusbar = self.setup_statusbar()
        self.dockwidget = self.setup_dockwidget()

        self.stacked_pages = self.setup_stacked_pages()
        self.toolbar = self.setup_toolbar()

        self.init_keyboard_shortcuts()
        self.connect_signals()

        MultiProcessWorker.create_pool()

    def connect_signals(self) -> None:
        """Connect signals to slots."""
        self.signals = Application.instance().signals

        self.signals.ask_download_database.connect(self.ask_download_database)

    def findChild(  # type: ignore [override] # explicitly override type hints
        self, type_: type["ChildType"], name: str, options: Optional["Qt.FindChildOption"] = None
    ) -> "ChildType":
        if options is None:
            options = Qt.FindChildOption.FindChildrenRecursively
        return super().findChild(type_, name, options)  # type: ignore [return-value] # explicitly override type hints

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

        dockwidget.setMinimumWidth(375)
        dockwidget.setVisible(False)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dockwidget)
        return dockwidget

    def setup_stacked_pages(self) -> NamedStackedWidget["BasePage"]:
        """Set up the different pages for each toolbar option."""
        stacked_pages = NamedStackedWidget["BasePage"]()
        self.setCentralWidget(stacked_pages)

        stacked_pages.addNamedWidget(OneAtomPage(), "OneAtomPage")
        stacked_pages.addNamedWidget(TwoAtomsPage(), "TwoAtomsPage")
        stacked_pages.addNamedWidget(LifetimesPage(), "LifetimesPage")
        # stacked_pages.addNamedWidget(C6Page(), "C6Page")

        # stacked_pages.addNamedWidget(SettingsPage(), "SettingsPage")
        # stacked_pages.addNamedWidget(AboutPage(), "AboutPage")
        return stacked_pages

    def setup_toolbar(self) -> QToolBar:
        """Set up the toolbar with icon buttons."""
        toolbar = QToolBar("Sidebar")
        toolbar.setMovable(False)
        toolbar.setOrientation(Qt.Orientation.Vertical)
        toolbar.setIconSize(QSize(32, 32))
        toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)

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

        default_page = "OneAtomPage"
        self.findChild(QAction, default_page).setChecked(True)
        self.stacked_pages.setCurrentNamedWidget(default_page)

        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, toolbar)

        return toolbar

    def init_keyboard_shortcuts(self) -> None:
        """Initialize keyboard shortcuts."""
        # Add Ctrl+W shortcut to close the window
        close_shortcut = QShortcut(QKeySequence("Ctrl+W"), self)
        close_shortcut.activated.connect(lambda: logger.info("Ctrl+W detected. Shutting down gracefully..."))
        close_shortcut.activated.connect(self.close)

    def closeEvent(self, event: QCloseEvent) -> None:
        """Make sure to also call Application.quit() when closing the window."""
        logger.debug("Close event triggered.")
        Application.quit()
        event.accept()

    def ask_download_database(self, species: str) -> bool:
        msg_box = QMessageBox()
        msg_box.setWindowTitle("Download missing database tables?")
        msg_box.setText(f"Database tables for {species} not found.")
        msg_box.setInformativeText("Would you like to download the missing database tables?")
        msg_box.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)

        download = msg_box.exec() == QMessageBox.StandardButton.Yes
        if download:
            self.statusbar.showMessage("Downloading database table ...", timeout=0)

            worker = MultiThreadWorker(lambda: download_databases_mp([species]))
            worker.enable_busy_indicator(self.stacked_pages.currentWidget())

            msg = "Successfully downloaded database table for " + species
            worker.signals.result.connect(lambda _result: self.statusbar.showMessage(msg, timeout=0))
            worker.signals.result.connect(lambda _result: setattr(Database, "_global_database", None))
            worker.signals.result.connect(lambda _result: MultiProcessWorker.terminate_all(create_new_pool=True))
            page = self.stacked_pages.currentWidget()
            if isinstance(page, SimulationPage):
                ket_config = page.ket_config
                for i in range(ket_config.n_atoms):
                    worker.signals.result.connect(
                        lambda _, atom=i: ket_config.signal_species_changed.emit(ket_config.get_species(atom), atom)
                    )
            worker.start()

        return download
