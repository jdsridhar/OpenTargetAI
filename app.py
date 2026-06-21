"""
OpenTargetAI - Open Source Ligand-Based Target Prediction Platform
Main application entry point.
"""

import sys
import logging
import os
from pathlib import Path

from PyQt6.QtWidgets import QApplication, QSplashScreen, QMessageBox
from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QPixmap, QColor, QPainter, QFont

# Ensure project root is in path
PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

from ui.main_window import MainWindow
from database.db_manager import DatabaseManager


def setup_logging() -> None:
    """Configure application-wide logging."""
    log_dir = PROJECT_ROOT / "logs"
    log_dir.mkdir(exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.FileHandler(log_dir / "opentargetai.log"),
            logging.StreamHandler(sys.stdout),
        ],
    )


def create_splash_screen() -> QSplashScreen:
    """Create a splash screen for app startup."""
    pixmap = QPixmap(600, 300)
    pixmap.fill(QColor("#0d1117"))

    painter = QPainter(pixmap)
    painter.setPen(QColor("#58a6ff"))

    title_font = QFont("Arial", 28, QFont.Weight.Bold)
    painter.setFont(title_font)
    painter.drawText(pixmap.rect(), Qt.AlignmentFlag.AlignCenter, "OpenTargetAI")

    sub_font = QFont("Arial", 12)
    painter.setFont(sub_font)
    painter.setPen(QColor("#8b949e"))
    rect = pixmap.rect()
    rect.setTop(rect.top() + 160)
    painter.drawText(rect, Qt.AlignmentFlag.AlignHCenter, "Open Source Ligand-Based Target Prediction")

    painter.end()

    splash = QSplashScreen(pixmap)
    splash.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint)
    return splash


def initialize_database() -> DatabaseManager:
    """Initialize and return the database manager."""
    db_path = PROJECT_ROOT / "data" / "opentargetai.db"
    db_path.parent.mkdir(exist_ok=True)
    db_manager = DatabaseManager(str(db_path))
    db_manager.initialize_schema()
    return db_manager


def main() -> None:
    """Main application entry point."""
    setup_logging()
    logger = logging.getLogger(__name__)
    logger.info("Starting OpenTargetAI...")

    app = QApplication(sys.argv)
    app.setApplicationName("OpenTargetAI")
    app.setApplicationVersion("1.0.0")
    app.setOrganizationName("OpenTargetAI")

    # Apply dark stylesheet
    app.setStyleSheet(load_stylesheet())

    # Show splash
    splash = create_splash_screen()
    splash.show()
    app.processEvents()

    try:
        # Initialize database
        splash.showMessage("Initializing database...", Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter, QColor("#58a6ff"))
        app.processEvents()
        db_manager = initialize_database()

        splash.showMessage("Loading application...", Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter, QColor("#58a6ff"))
        app.processEvents()

        window = MainWindow(db_manager)
        window.show()

        # Close the splash as soon as the main window is up, then (once the
        # event loop is running and nothing is on top) offer the demo dataset.
        splash.finish(window)
        QTimer.singleShot(300, window.offer_demo_if_empty)

        logger.info("Application started successfully.")
        sys.exit(app.exec())

    except Exception as e:
        logger.exception("Fatal error during startup")
        splash.close()
        QMessageBox.critical(None, "Startup Error", f"Failed to start OpenTargetAI:\n{e}")
        sys.exit(1)


def load_stylesheet() -> str:
    """Load the dark scientific theme stylesheet."""
    return """
    QMainWindow, QDialog, QWidget {
        background-color: #0d1117;
        color: #c9d1d9;
        font-family: 'Segoe UI', Arial, sans-serif;
        font-size: 13px;
    }
    QTabWidget::pane {
        border: 1px solid #30363d;
        background-color: #161b22;
    }
    QTabBar::tab {
        background-color: #21262d;
        color: #8b949e;
        padding: 8px 16px;
        border: 1px solid #30363d;
        border-bottom: none;
        margin-right: 2px;
    }
    QTabBar::tab:selected {
        background-color: #161b22;
        color: #58a6ff;
        border-top: 2px solid #58a6ff;
    }
    QTabBar::tab:hover {
        background-color: #30363d;
        color: #c9d1d9;
    }
    QPushButton {
        background-color: #21262d;
        color: #c9d1d9;
        border: 1px solid #30363d;
        border-radius: 6px;
        padding: 6px 14px;
        font-weight: bold;
    }
    QPushButton:hover {
        background-color: #30363d;
        border-color: #58a6ff;
        color: #58a6ff;
    }
    QPushButton:pressed {
        background-color: #1f6feb;
        color: #ffffff;
    }
    QPushButton#runButton {
        background-color: #1f6feb;
        color: #ffffff;
        font-size: 14px;
        padding: 10px 20px;
        border-radius: 8px;
    }
    QPushButton#runButton:hover {
        background-color: #388bfd;
    }
    QPushButton#runButton:disabled {
        background-color: #21262d;
        color: #484f58;
    }
    QLineEdit, QTextEdit, QPlainTextEdit {
        background-color: #21262d;
        color: #c9d1d9;
        border: 1px solid #30363d;
        border-radius: 6px;
        padding: 6px;
        selection-background-color: #1f6feb;
    }
    QLineEdit:focus, QTextEdit:focus {
        border-color: #58a6ff;
    }
    QComboBox {
        background-color: #21262d;
        color: #c9d1d9;
        border: 1px solid #30363d;
        border-radius: 6px;
        padding: 6px 10px;
    }
    QComboBox::drop-down {
        border: none;
    }
    QComboBox QAbstractItemView {
        background-color: #21262d;
        color: #c9d1d9;
        border: 1px solid #30363d;
        selection-background-color: #1f6feb;
    }
    QSlider::groove:horizontal {
        height: 4px;
        background: #30363d;
        border-radius: 2px;
    }
    QSlider::handle:horizontal {
        background: #58a6ff;
        width: 14px;
        height: 14px;
        margin: -5px 0;
        border-radius: 7px;
    }
    QSlider::sub-page:horizontal {
        background: #1f6feb;
        border-radius: 2px;
    }
    QTableWidget {
        background-color: #161b22;
        color: #c9d1d9;
        border: 1px solid #30363d;
        gridline-color: #21262d;
        alternate-background-color: #1c2128;
    }
    QTableWidget::item:selected {
        background-color: #1f6feb;
        color: #ffffff;
    }
    QHeaderView::section {
        background-color: #21262d;
        color: #8b949e;
        padding: 8px;
        border: none;
        border-bottom: 1px solid #30363d;
        font-weight: bold;
    }
    QHeaderView::section:hover {
        background-color: #30363d;
        color: #c9d1d9;
    }
    QScrollBar:vertical {
        background: #0d1117;
        width: 10px;
    }
    QScrollBar::handle:vertical {
        background: #30363d;
        border-radius: 5px;
        min-height: 20px;
    }
    QScrollBar::handle:vertical:hover {
        background: #58a6ff;
    }
    QScrollBar:horizontal {
        background: #0d1117;
        height: 10px;
    }
    QScrollBar::handle:horizontal {
        background: #30363d;
        border-radius: 5px;
        min-width: 20px;
    }
    QProgressBar {
        background-color: #21262d;
        border: 1px solid #30363d;
        border-radius: 4px;
        text-align: center;
        color: #c9d1d9;
    }
    QProgressBar::chunk {
        background-color: #1f6feb;
        border-radius: 4px;
    }
    QGroupBox {
        border: 1px solid #30363d;
        border-radius: 8px;
        margin-top: 12px;
        padding-top: 8px;
        color: #8b949e;
        font-weight: bold;
    }
    QGroupBox::title {
        subcontrol-origin: margin;
        subcontrol-position: top left;
        padding: 0 6px;
        color: #58a6ff;
    }
    QLabel#sectionTitle {
        color: #58a6ff;
        font-size: 14px;
        font-weight: bold;
    }
    QLabel#statValue {
        color: #3fb950;
        font-size: 13px;
        font-weight: bold;
    }
    QSplitter::handle {
        background: #30363d;
        border-radius: 3px;
    }
    QSplitter::handle:horizontal {
        width: 8px;
        margin: 3px 1px;
    }
    QSplitter::handle:vertical {
        height: 8px;
        margin: 1px 3px;
    }
    QSplitter::handle:hover {
        background: #58a6ff;
    }
    QSplitter::handle:pressed {
        background: #1f6feb;
    }
    QStatusBar {
        background-color: #161b22;
        color: #8b949e;
        border-top: 1px solid #30363d;
    }
    QMenuBar {
        background-color: #161b22;
        color: #c9d1d9;
        border-bottom: 1px solid #30363d;
    }
    QMenuBar::item:selected {
        background-color: #30363d;
    }
    QMenu {
        background-color: #21262d;
        color: #c9d1d9;
        border: 1px solid #30363d;
    }
    QMenu::item:selected {
        background-color: #1f6feb;
    }
    QToolTip {
        background-color: #21262d;
        color: #c9d1d9;
        border: 1px solid #30363d;
        border-radius: 4px;
        padding: 4px;
    }
    QSpinBox, QDoubleSpinBox {
        background-color: #21262d;
        color: #c9d1d9;
        border: 1px solid #30363d;
        border-radius: 6px;
        padding: 4px;
    }
    QCheckBox {
        color: #c9d1d9;
        spacing: 6px;
    }
    QCheckBox::indicator {
        width: 16px;
        height: 16px;
        border: 1px solid #30363d;
        border-radius: 3px;
        background-color: #21262d;
    }
    QCheckBox::indicator:checked {
        background-color: #1f6feb;
        border-color: #58a6ff;
    }
    QRadioButton {
        color: #c9d1d9;
        spacing: 6px;
    }
    QRadioButton::indicator {
        width: 14px;
        height: 14px;
        border: 1px solid #30363d;
        border-radius: 7px;
        background-color: #21262d;
    }
    QRadioButton::indicator:checked {
        background-color: #1f6feb;
        border-color: #58a6ff;
    }
    """


if __name__ == "__main__":
    main()
