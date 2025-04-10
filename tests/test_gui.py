import time
from typing import TYPE_CHECKING

from pairinteraction_gui.main_window import MainWindow

if TYPE_CHECKING:
    from pytestqt.qtbot import QtBot

    from pairinteraction_gui.app import Application
    from pairinteraction_gui.page import OneAtomPage


def test_main_window_basic(qtbot: "QtBot", qapp: "Application") -> None:
    """Test basic main window functionality."""
    window = MainWindow()
    window.show()
    qtbot.addWidget(window)

    one_atom_page: OneAtomPage = window.stacked_pages.getNamedWidget("OneAtomPage")
    ket_label = one_atom_page.ket_config.ket_label[0].text()
    assert all(x in ket_label for x in ["Rb", "60", "S", "1/2"])

    qn_item = one_atom_page.ket_config.stacked_qn[0].currentWidget().items[0]
    assert qn_item.label == "n"
    assert qn_item.value() == 60
    qn_item.spinboxes["value"].setValue(61)
    assert qn_item.value() == 61
    ket_label = one_atom_page.ket_config.ket_label[0].text()
    assert all(x in ket_label for x in ["Rb", "61", "S", "1/2"])

    # one_atom_page.calculate_and_abort.getNamedWidget("Calculate").click()
    # time.sleep(1)  # Simulate some processing time
    # one_atom_page.calculate_and_abort.getNamedWidget("Abort").click()
    # time.sleep(2)  # Simulate some processing time
