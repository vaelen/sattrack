#include <sattrack/gui.h>

class SatTrackWindow : public Gtk::Window {
public:
    SatTrackWindow() {
        set_title("SatTrack - Satellite Tracker");
        set_default_size(400, 300);

        // Add your UI components here
        m_label.set_text("SatTrack - Satellite tracking application");
        add(m_label);

        show_all_children();
    }

    virtual ~SatTrackWindow() {}

protected:
    Gtk::Label m_label;
};

int startGUI(int argc, char* argv[]) {
    auto app = Gtk::Application::create(argc, argv, "org.sattrack.app");

    SatTrackWindow window;

    return app->run(window);
}
