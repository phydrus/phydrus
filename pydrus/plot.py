import matplotlib.pyplot as plt
from matplotlib import cm


class Plots:
    def __init__(self, ml):
        self.ml = ml

    def profile(self, figsize=(4, 10), title="Soil Profile", cmap="YlOrBr",
                color_by="Ks", **kwargs):
        """Method to plot the soil profile.

        Parameters
        ----------
        figsize: tuple, optional
        title: str, optional
        cmap: str, optional
            String with a named Matplotlib colormap.
        color_by: str, optional
            Column from the material properties sed to color the materials.
            Default is "Ks".

        Returns
        -------
        ax: matplotlib axes instance

        """
        fig, ax = plt.subplots(figsize=figsize, **kwargs)

        top = self.ml.profile.loc[:, "x"].max()
        w = self.ml.profile.loc[:, "h"].max()
        w = w + 0.2 * w

        # Set colors by color_by
        col = self.ml.materials[color_by]
        col = (col - col.min()) / (col.max() - col.min())
        colors = cm.get_cmap(cmap, 7)(col.values)

        for i in self.ml.profile.index[1:]:
            bot = self.ml.profile.loc[i, "x"]
            h = bot - top
            color = colors[self.ml.profile.loc[i, "Mat"] - 1]
            patch = plt.Rectangle(xy=(0, top), width=w, height=h, linewidth=1,
                                  edgecolor="darkgray", facecolor=color)
            ax.add_patch(patch)
            top = bot

        line = ax.plot(self.ml.profile.loc[:, ["h"]].values,
                       self.ml.profile.loc[:, ["x"]].values,
                       label="Initial head")

        ax.set_xlim(0, w)
        ax.set_ylim(self.ml.profile.loc[:, "x"].min(),
                    self.ml.profile.loc[:, "x"].max())
        ax.set_xlabel("h [{}]".format(self.ml.basic_information["LUnit"]))
        ax.set_ylabel("depth [{}]".format(self.ml.basic_information["LUnit"]))
        ax.set_title(title)

        legend_elements = [line[0]]
        for i, color in enumerate(colors):
            legend_elements.append(plt.Rectangle((0, 0), 0, 0, color=color,
                                                 label="material {}".format(
                                                     i)))

        plt.legend(handles=legend_elements, loc="best")

        plt.tight_layout()
        return ax
    
    def profile_information(self, figsize=(6, 6), title="Profile Information",
                            cmap="YlOrBr", **kwargs):
        """Method to plot the soil profile.

        Parameters
        ----------
        figsize: tuple, optional
        title: str, optional
        cmap: str, optional
            String with a named Matplotlib colormap.
        color_by: str, optional
            Column from the material properties sed to color the materials.
            Default is "Ks".

        Returns
        -------
        ax: matplotlib axes instance

        """        
        df = self.ml.read_nod_inf()[self.ml.time_information["tMax"]]
        df1 = self.ml.read_nod_inf()[self.ml.time_information["tInit"]]
        l_unit = self.ml.basic_information["LUnit"]
        t_unit = self.ml.basic_information["TUnit"]
        use_cols = ("Head","Moisture","K", "C", "Flux", "Sink")
        col_names = ("Pressure Head", "Water Content", "Hydraulic Conductivity", 
                 "Hydraulic Capacity", "Water Flux", "Root Uptake")
        units = ["h [{}]".format(l_unit),
                  "Theta [-]","K [{}/days]".format(l_unit),
                  "C [1/{}]".format(l_unit), 
                  "v [{}/{}]".format(l_unit, t_unit),
                  "S [1/{}]".format(t_unit)]
        
        for cols,names,unit in zip(use_cols, col_names,units):
            fig, ax = plt.subplots(figsize=figsize, **kwargs)
            if cols == "Moisture":
                plt.axvline(df1[cols][1],color = 'k')                          
            if cols == "K":
                plt.axvline(df1[cols][1],color='k')
            if cols == "C":
                plt.axvline(df1[cols][1],color='k')
            if cols == "Flux" and t_unit == "min":
                df[cols] = df[cols]*60
                unit = "v [{}/days]".format(l_unit)

            line = ax.plot(df[cols], self.ml.profile.loc[:, ["x"]].values, "b")   
            ax.set_ylim(self.ml.profile.loc[:, "x"].min(),self.ml.profile.loc[:, "x"].max())
            ax.set_xlabel(unit)
            ax.set_ylabel("Depth [{}]".format(self.ml.basic_information["LUnit"]))
            ax.set_title("Profile Information: " + names)
            ax.grid(linestyle='--')
            plt.tight_layout()
        return ax
    def mass_balance(self, figsize=(6, 10), title="Mass_Balance_Information",
                        **kwargs):
        """Method to show the Soil Hydraulic Properties.

        Parameters
        ----------
        figsize: tuple, optional
        title: str, optional

        Returns
        -------
        balance: opens BALANCE.OUT in notepad
        """
        folder = self.ml.ws_name
        import sys
        import subprocess
        balance = subprocess.Popen(["notepad.exe",folder + "\\BALANCE.OUT"])
        return balance