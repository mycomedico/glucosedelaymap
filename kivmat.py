# importing pyplot for graph plotting
from matplotlib import pyplot as plt

# importing numpy
import numpy as np
from backend_kivyagg import FigureCanvasKivyAgg

# importing kivyapp
from kivy.app import App

# importing kivy builder
from kivy.lang import Builder


# this is the main class which will
# render the whole application
class uiApp(App):

    def build(self):
        self.str = Builder.load_string(""" 

BoxLayout:
    layout:layout

    BoxLayout:

        id:layout

                                """)

        signal = [7, 89.6, 45. - 56.34]

        signal = np.array(signal)

        # this will plot the signal on graph
        plt.plot(signal)

        # setting x label
        plt.xlabel('Time(s)')

        # setting y label
        plt.ylabel('signal (norm)')
        plt.grid(True, color='lightgray')

        # adding plot to kivy boxlayout
        self.str.layout.add_widget(FigureCanvasKivyAgg(plt.gcf()))
        return self.str


# running the application
uiApp().run()