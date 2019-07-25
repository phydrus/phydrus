import matplotlib.pyplot as plt


class Plots:
    def __init__(self, ml):
        self.ml = ml

    def profile(self, figsize=(3, 10), title="Soil Profile", **kwargs):
        fig, ax = plt.subplots(figsize=figsize, **kwargs)

        top = self.ml.profile.loc[:, "x"].max()
        w  = self.ml.profile.loc[:, "h"].max()
        w = w + 0.1 * w
        c = ["C0", "C1", "C2"]

        for i in self.ml.profile.index[1:]:
            bot = self.ml.profile.loc[i, "x"]
            h = bot - top
            color = c[self.ml.profile.loc[i, "Mat"]]
            patch = plt.Rectangle(xy=(0, top), width=w, height=h, linewidth=1,
                                  edgecolor="k", facecolor=color)
            ax.add_patch(patch)
            top = bot

        ax.plot(self.ml.profile.loc[:, ["h"]].values,
                self.ml.profile.loc[:, ["x"]].values)

        ax.set_xlim(0, w)
        ax.set_ylim(self.ml.profile.loc[:, "x"].min(),
                    self.ml.profile.loc[:, "x"].max())
        ax.set_xlabel("h [{}]".format(self.ml.basic_information["LUnit"]))
        ax.set_ylabel("depth [{}]".format(self.ml.basic_information["LUnit"]))
        ax.set_title(title)
        plt.tight_layout()
        return ax

