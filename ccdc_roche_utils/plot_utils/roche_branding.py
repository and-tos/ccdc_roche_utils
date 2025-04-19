from matplotlib.colors import LinearSegmentedColormap
import numpy as np


class RocheColours(object):
    def __init__(self):
        roche_colours_rgb_dict = {
            "roche_blue": (11, 65, 205),
            "dark_blue": (2, 35, 102),
            "light_blue": (20, 130, 250),
            "extra_light_blue": (189, 227, 255),
            "neutral1": (250, 201, 181),
            "neutral2": (250, 214, 199),
            "neutral3": (255, 232, 222),
            "neutral4": (255, 247, 245),
            "extra_dark_red": (140, 0, 0),
            "dark_red": (196, 0, 0),
            "red": (255, 31, 38),
            "light_red": (255, 135, 130),
            "extra_dark_orange": (178, 43, 13),
            "dark_orange": (237, 74, 13),
            "orange": (255, 125, 41),
            "light_orange": (255, 189, 105),
            "extra_dark_purple": (125, 0, 150),
            "dark_purple": (188, 54, 240),
            "purple": (224, 133, 252),
            "light_purple": (242, 212, 255),
        }
        for key in roche_colours_rgb_dict.keys():
            roche_colours_rgb_dict[key] = tuple(
                i / 255 for i in roche_colours_rgb_dict[key]
            )
        self.roche_colours_rgb_dict = roche_colours_rgb_dict
        self.cm = LinearSegmentedColormap.from_list(
            "roche",
            list(roche_colours_rgb_dict.values()),
            N=len(roche_colours_rgb_dict),
        )
        self.palette = self.cm(np.linspace(0, 1, len(roche_colours_rgb_dict)))
