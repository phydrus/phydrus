import matplotlib.pyplot as plt
from matplotlib import cm


class Plots:
    def __init__(self, ml):
        self.ml = ml

    def profile(self, figsize=(3, 6), title="Soil Profile", cmap="YlOrBr",
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
        col = self.ml.materials["water"][color_by]
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
        ax.set_xlabel("h [{}]".format(self.ml.basic_info["LUnit"]))
        ax.set_ylabel("depth [{}]".format(self.ml.basic_info["LUnit"]))
        ax.set_title(title)

        legend_elements = [line[0]]
        for i, color in enumerate(colors):
            legend_elements.append(plt.Rectangle((0, 0), 0, 0, color=color,
                                                 label="material {}".format(i))
                                   )

        plt.legend(handles=legend_elements, loc="best")
        plt.tight_layout()
        return ax

    def profile_information(self, data="Pressure Head", times=None,
                            legend=True, figsize=(5, 3), **kwargs):
        """Method to plot the soil profile information.

        Parameters
        ----------
        data: str, optional
            String with the variable of the profile information to plot. 
            You can choose between: "Pressure Head", "Water Content", 
            "Hydraulic Conductivity","Hydraulic Capacity", "Water Flux", 
            "Root Uptake".
            Default is "Pressure Head".
        times: list of int
            List of integers of the time step to plot.       
        figsize: tuple, optional
        legend: boolean, optional

        Returns
        -------
        ax: matplotlib axes instance

        """
        l_unit = self.ml.basic_info["LUnit"]
        t_unit = self.ml.basic_info["TUnit"]
        m_unit = self.ml.basic_info["MUnit"]

        use_cols = ("Head", "Moisture", "K", "C", "Flux", "Sink")
        col_names = ("Pressure Head", "Water Content",
                     "Hydraulic Conductivity", "Hydraulic Capacity",
                     "Water Flux", "Root Uptake")
        units = ["h [{}]".format(l_unit),
                 "Theta [-]", "K [{}/days]".format(l_unit),
                 "C [1/{}]".format(l_unit),
                 "v [{}/{}]".format(l_unit, t_unit),
                 "S [1/{}]".format(t_unit)]

        if self.ml.basic_info["lChem"]:
            use_cols = use_cols + ("Conc(1..NS)", "Sorb(1...NS)")
            col_names = col_names + ("Concentration", "Sorbtion")
            units.extend(["c [{}/{}*3]".format(m_unit, l_unit), "sorb."])

        col = col_names.index(data)
        fig, ax = plt.subplots(figsize=figsize, **kwargs)
        dfs = self.ml.read_nod_inf(times=times)

        if times is None or len(times) > 1:
            for key, df in dfs.items():
                df.plot(x=use_cols[col], y="Depth", ax=ax,
                        label="time= " + str(key))
        else:
            dfs.plot(x=use_cols[col], y="Depth", ax=ax,
                     label="T " + str(times))

        ax.set_xlabel(units[col])
        ax.set_ylabel("Depth [{}]".format(self.ml.basic_info["LUnit"]))
        ax.grid(linestyle='--')

        if legend:
            ax.legend(bbox_to_anchor=(1, 1), loc="upper left")

        plt.tight_layout()
        return ax

    def water_flow(self, data="Potential Surface Flux", figsize=(6, 3),
                   **kwargs):
        """Method to plot the water flow information.

        Parameters
        ----------
        data: str, optional
            String with the variable of the water flow information to plot.
            You can choose between: "Potential Surface Flux",
            "Potential Root Water Uptake", "Actual Surface Flux", 
            "Actual Root Water Uptake", "Bottom Flux", 
            "Pressure head at the soil surface", 
            "Mean value of the pressure head over the region",
            "Pressure head at the Bottom of the soil profile",
            "Surface runoff", "Volume of water in the entire flow domain".
            Default is "Potential Surface Flux".
        figsize: tuple, optional

        Returns
        -------
        ax: matplotlib axes instance

        """
        col_names = ("Potential Surface Flux", "Potential Root Water Uptake",
                     "Actual Surface Flux", "Actual Root Water Uptake",
                     "Bottom Flux", "Pressure head at the soil surface",
                     "Mean value of the pressure head over the region",
                     "Pressure head at the Bottom of the soil profile",
                     "Surface runoff",
                     "Volume of water in the entire flow domain")

        cols = ("rTop", "rRoot", "vTop", "vRoot", "vBot", "hTop", "hRoot",
                "hBot", "RunOff", "Volume")
        df = self.ml.read_tlevel()
        col = col_names.index(data)

        if col < 5:
            fig, axes = plt.subplots(1, 2, figsize=figsize, **kwargs)
            df.plot(y=cols[col], ax=axes[0], use_index=True)
            axes[0].set_ylabel(data)
            axes[0].set_xlabel("Time [{}]".format(self.ml.basic_info["TUnit"]))

            # Cumulative sum
            df.plot(y="sum(" + cols[col] + ")", ax=axes[1], use_index=True)
            axes[1].set_ylabel("Cum. {}".format(data))
            axes[1].set_xlabel("Time [{}]".format(self.ml.basic_info["TUnit"]))
            fig.tight_layout()
        else:
            fig, axes = plt.subplots(1, 1, figsize=figsize, **kwargs)
            axes.plot(df.index, df[cols[col]])
            axes.set_ylabel(data)
            axes.set_xlabel("Time [{}]".format(self.ml.basic_info["TUnit"]))
        return axes

    def soil_properties(self, data="Water Content", figsize=(6, 3), **kwargs):
        """Method to plot the soil hydraulic properties.

        Parameters
        ----------
        data: str, optional
            String with the variable of the water flow information to plot. 
            You can choose between: "Water Content", "Pressure head",
            "log Pressure head", "Hydraulic Capacity", "Hydraulic
            Conductivity", "log Hydraulic Conductivity", "Effective Water
            Content". Default is "Water Content".
        figsize: tuple, optional

        Returns
        -------
        axes: matplotlib axes instance

        """
        col_names = ("Water Content", "Pressure head", "log Pressure head",
                     "Hydraulic Capacity", "Hydraulic Conductivity",
                     "log Hydraulic Conductivity", "Effective Water Content")
        cols = ("theta", "h", "log_h", "C", "K", "log_K", "S", "Kv")
        col = col_names.index(data)

        dfs = self.ml.read_i_check()

        _, axes = plt.subplots(figsize=figsize, nrows=1, ncols=2,
                               sharey=True, **kwargs)

        for i, df in dfs.items():
            name = "Node {}".format(i)
            df.plot(x="h", y=cols[col], ax=axes[0], label=name)
            df.plot(x="log_h", y=cols[col], ax=axes[1], label=name)

        axes[0].set_xlabel(xlabel="h")
        axes[1].set_xlabel(xlabel="log_h")
        axes[0].set_ylabel(cols[col])

        return axes

    def obs_points(self, data="Pressure Head", figsize=(4, 3), **kwargs):
        """Method to plot the pressure heads, water contents and water fluxes
        at selected observation nodes.

        Parameters
        ----------
        data: str, optional
            String with the variable of the variable to plot.
            You can choose between: "Pressure head", "Water content",
            "Water flux".
        figsize: tuple, optional

        Returns
        -------
        axes: matplotlib axes instance

        """
        col_names = ("Pressure Head", "Water Content", "Water Flux",
                     "Concentration")
        cols = ("h", "theta", "Temp", "Conc")
        col = col_names.index(data)

        dfs = self.ml.read_obs_node()

        _, ax = plt.subplots(figsize=figsize, **kwargs)
        for i, df in dfs.items():
            name = "Node {}".format(i)
            df.plot(y=cols[col], ax=ax, label=name, use_index=True)

        ax.set_xlabel("Time [{}]".format(self.ml.basic_info["TUnit"]))
        ax.set_ylabel(cols[col])
        return ax
