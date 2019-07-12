import matplotlib.pyplot as plt


class Plots:
    def __init__(self, ml):
        self.ml = ml

    def plot_profile(self, figsize=(3, 10), **kwargs):
        fig, axes = plt.subplots(figsize, **kwargs)



        return axes
