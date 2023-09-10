import sys
from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtWidgets import QApplication, QMainWindow, QToolBar, QLineEdit, QPushButton, QVBoxLayout, QWidget, QSizePolicy, QProgressBar
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEnginePage

class BrowserWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.browser = QWebEngineView()
        self.browser.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.browser.setUrl(QUrl("https://duckduckgo.com"))  # Set DuckDuckGo as the initial URL
        self.setCentralWidget(self.browser)

        # Set dark theme (you can customize the colors)
        dark_theme = "QWebEngineView { background-color: #121212; color: #FFFFFF; }"
        self.browser.setStyleSheet(dark_theme)

        # Create a navigation toolbar
        nav_toolbar = QToolBar("Navigation")
        self.addToolBar(nav_toolbar)

        # Add an address bar for entering URLs or search queries
        self.address_bar = QLineEdit()
        nav_toolbar.addWidget(self.address_bar)

        # Add a "Go" button to navigate to the entered URL or perform searches
        go_button = QPushButton("Go/Search")
        go_button.clicked.connect(self.navigate_or_search)
        nav_toolbar.addWidget(go_button)

        # Connect Enter key press to navigate or search
        self.address_bar.returnPressed.connect(self.navigate_or_search)

        # Add a progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.progress_bar.setVisible(False)
        nav_toolbar.addWidget(self.progress_bar)

        # Add a reload button
        reload_button = QPushButton("Reload")
        reload_button.clicked.connect(self.reload_page)
        nav_toolbar.addWidget(reload_button)

        # Add an abort button
        abort_button = QPushButton("Abort")
        abort_button.clicked.connect(self.abort_loading)
        nav_toolbar.addWidget(abort_button)

        # Add forward and backward buttons
        nav_toolbar.addSeparator()
        back_button = QPushButton("Back")
        back_button.clicked.connect(self.browser.back)
        nav_toolbar.addWidget(back_button)

        forward_button = QPushButton("Forward")
        forward_button.clicked.connect(self.browser.forward)
        nav_toolbar.addWidget(forward_button)

        # Connect the progress bar to page load events
        self.browser.page().loadStarted.connect(self.start_loading)
        self.browser.page().loadFinished.connect(self.finish_loading)

        # Set the browser window title
        self.setWindowTitle("Humsafar")

    def navigate_or_search(self):
        query = self.address_bar.text()
        if not query.startswith(("http://", "https://")):
            # If it's not a URL, perform a DuckDuckGo search
            search_url = f"https://duckduckgo.com/?q={query.replace(' ', '+')}"
            self.browser.setUrl(QUrl(search_url))
        else:
            # If it's a URL, navigate to the entered URL
            self.browser.setUrl(QUrl(query))

    def start_loading(self):
        # Show the progress bar when loading starts
        self.progress_bar.setVisible(True)

    def finish_loading(self):
        # Hide the progress bar when loading is finished
        self.progress_bar.setVisible(False)

    def reload_page(self):
        # Reload the current page
        self.browser.reload()

    def abort_loading(self):
        # Abort the current page loading
        self.browser.page().triggerAction(QWebEnginePage.Stop)

def main():
    app = QApplication(sys.argv)
    window = BrowserWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
