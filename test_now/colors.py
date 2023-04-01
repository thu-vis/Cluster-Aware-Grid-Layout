
class MyColorMap(object):
    def __init__(self) -> None:
        self.rgbs = [
            (171, 227, 246),
            (248, 182, 187),
            (185, 243, 190),
            (243, 228, 191),
            (216, 193, 246),
            (252, 245, 155),
            (221, 221, 221),
            (138, 170, 208),
            (191, 185, 134),
            (255, 193, 152),
            (127, 137, 253),
            (255, 136, 104),
            (175, 203, 191),
            (170, 167, 188),
            (254, 228, 179)
        ]
        self.rgbs1 = list(map(lambda x: list(map(lambda v: v / 255, x)), self.rgbs))
    
    def colorSet(self, num):
        return self.rgbs1[:num]

    def color(self, num):
        if num < 0:
            return (1, 1, 1)
        return self.rgbs1[num]
