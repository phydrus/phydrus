import matplotlib.pyplot as plt


class Plots:
    """Class that contains all the methods to plot a Phydrus Model.

    Parameters
    ----------
    ml: phydrus.Model
        Phydrus Model Instance to connect the methods to the model.

    Examples
    --------
    ml.plots.profile()

    """

    def __init__(self, ml):
        self.ml = ml

    def profile(
        self,
        figsize=(3, 6),
        cmap="YlOrBr_r",
        color_by="Ks",
        vmin=None,
        vmax=None,
        show_grid=True,
        plot_h=True,
        ax=None,
    ):
        """Method to plot the soil profyle

        Parameters
        ----------
        figsize : tuple, optional
            Tuple with the size of the figure, by default (3, 6)
        cmap : str, optional
            String with a named Matplotlib colormap, by default 'YlOrBr_r'
        color_by : str, optional
            Column from the material properties to color by, by default 'Ks'
        vmin : str, optional
            Minimum value for the colormap, if None use minimum of color_by
        vmax : str, optional
            Maximum value for the colormap, if None use minimum of color_by
        show_grid : bool, optional
            Show the grid lines of dx, by default True
        ax : Matplotlib.Axes, optional
            Axes to plot the profile in, by default None which creates an Axes

        Returns
        -------
        ax: matplotlib axes instance

        """

        if ax is None:
            _, ax = plt.subplots(figsize=figsize)

        y = self.ml.profile.loc[:, "x"]
        top = y.max()
        hmax = self.ml.profile.loc[:, "h"].max()
        hmin = self.ml.profile.loc[:, "h"].min()

        # Set colors by color_by
        col = self.ml.materials["water"][color_by]
        if vmin is None:
            vmin = col.min()
        if vmax is None:
            vmax = col.max()
        coln = (col - vmin) / (vmax - vmin)
        colors = plt.cm.get_cmap(cmap, len(col))(coln.values)

        if show_grid:
            edgecolor = (0.169, 0.169, 0.169, 0.2)
        else:
            edgecolor = None

        for i in self.ml.profile.index[1:]:
            bot = y.loc[i]
            h = bot - top
            j = self.ml.profile.loc[i, "Mat"] - 1
            patch = plt.Rectangle(
                xy=(-1e6, top),
                width=1e11,
                height=h,
                linewidth=1,
                edgecolor=edgecolor,
                facecolor=colors[j],
                label=f"Material {j}",
            )
            ax.add_patch(patch)
            top = bot

        if plot_h:
            ax.plot(self.ml.profile.loc[:, ["h"]].values, y.values, label=r"$\psi_i$")

        ax.set_xlim(hmin - 0.2 * abs(hmin), hmax + 0.2 * abs(hmax))
        ax.set_ylim(y.min(), y.max())

        return ax

    def profile_information(
        self,
        data="Pressure Head",
        times=None,
        legend=True,
        figsize=(5, 3),
        ax=None,
        **kwargs,
    ):
        """
        Method to plot the soil profile information.

        Parameters
        ----------
        data: str, optional
            String with the variable of the profile information to plot.
            You can choose between: "Pressure Head", "Water Content",
            "Hydraulic Conductivity","Hydraulic Capacity", "Water Flux",
            "Root Uptake", "Temperature". Default is "Pressure Head".
        times: list of int
            List of integers of the time step to plot.
        figsize: tuple, optional
        legend: boolean, optional
        ax : Matplotlib.Axes, optional
            Axes to plot the profile in, by default None which creates an Axes

        Returns
        -------
        ax: matplotlib axes instance

        """
        l_unit = self.ml.basic_info["LUnit"]
        t_unit = self.ml.basic_info["TUnit"]
        m_unit = self.ml.basic_info["MUnit"]

        use_cols = ("Head", "Moisture", "K", "C", "Flux", "Sink", "Temp")

        col_names = (
            "Pressure Head",
            "Water Content",
            "Hydraulic Conductivity",
            "Hydraulic Capacity",
            "Water Flux",
            "Root Uptake",
            "Temperature",
        )
        units = [
            f"h [{l_unit}]",
            "Theta [-]",
            f"K [{l_unit}/days]",
            f"C [1/{l_unit}]",
            f"v [{l_unit}/{t_unit}]",
            f"S [1/{t_unit}]",
            "T [°C]",
        ]

        if self.ml.basic_info["lChem"]:
            use_cols = use_cols + ("Conc(1..NS)", "Sorb(1...NS)")
            col_names = col_names + ("Concentration", "Sorbtion")
            units.extend([f"c [{m_unit}/{l_unit}*3]", "sorb."])

        col = col_names.index(data)
        if ax is None:
            _, ax = plt.subplots(figsize=figsize, **kwargs)
        dfs = self.ml.read_nod_inf(times=times)

        if times is None or len(times) > 1:
            for key, df in dfs.items():
                df.plot(x=use_cols[col], y="Depth", ax=ax, label=f"time={key}")
        else:
            dfs.plot(x=use_cols[col], y="Depth", ax=ax, label=f"T {times}")

        ax.set_xlabel(units[col])
        ax.set_ylabel(f"Depth [{self.ml.basic_info['LUnit']}]")
        ax.grid(linestyle="--")

        if legend:
            ax.legend(bbox_to_anchor=(1, 1), loc="upper left")

        plt.tight_layout()
        return ax

    def water_flow(
        self, data="Potential Surface Flux", figsize=(6, 3), ax=None, **kwargs
    ):

        """
        Method to plot the water flow information.

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
        col_names = (
            "Potential Surface Flux",
            "Potential Root Water Uptake",
            "Actual Surface Flux",
            "Actual Root Water Uptake",
            "Bottom Flux",
            "Pressure head at the soil surface",
            "Mean value of the pressure head over the region",
            "Pressure head at the Bottom of the soil profile",
            "Surface runoff",
            "Volume of water in the entire flow domain",
        )

        cols = (
            "rTop",
            "rRoot",
            "vTop",
            "vRoot",
            "vBot",
            "hTop",
            "hRoot",
            "hBot",
            "RunOff",
            "Volume",
        )
        df = self.ml.read_tlevel()
        col = col_names.index(data)

        if ax is None:
            _, ax = plt.subplots(1, 2, figsize=figsize, **kwargs)
        df.plot(y=cols[col], ax=ax[0], use_index=True)
        ax[0].set_ylabel(data)
        ax[0].set_xlabel(f"Time [{self.ml.basic_info['TUnit']}]")

        # Cumulative sum
        df.plot(y=f"sum({cols[col]})", ax=ax[1], use_index=True)
        ax[1].set_ylabel(f"Cum. {data}")
        ax[1].set_xlabel(f"Time [{self.ml.basic_info['TUnit']}]")

        plt.tight_layout()
        return ax

    def soil_properties(self, data="Water Content", figsize=(6, 3), ax=None, **kwargs):
        """
        Method to plot the soil hydraulic properties.

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
        col_names = (
            "Water Content",
            "Pressure head",
            "log Pressure head",
            "Hydraulic Capacity",
            "Hydraulic Conductivity",
            "log Hydraulic Conductivity",
            "Effective Water Content",
        )
        cols = ("theta", "h", "log_h", "C", "K", "log_K", "S", "Kv")
        col = col_names.index(data)

        dfs = self.ml.read_i_check()

        if ax is None:
            _, ax = plt.subplots(
                figsize=figsize, nrows=1, ncols=2, sharey=True, **kwargs
            )
        for i, df in dfs.items():
            name = f"Node {i}"
            df.plot(x="h", y=cols[col], ax=ax[0], label=name)
            df.plot(x="log_h", y=cols[col], ax=ax[1], label=name)

        ax[0].set_xlabel(xlabel="h")
        ax[1].set_xlabel(xlabel="log_h")
        ax[0].set_ylabel(cols[col])

        return ax

    def obs_points(self, data="h", figsize=(4, 3), ax=None, **kwargs):
        """
        Method to plot the pressure heads, water contents and water fluxes.

        Parameters
        ----------
        data: str, optional
            String with the variable of the variable to plot.
            You can choose between: "h", "theta", "Temp", "Conc".
        figsize: tuple, optional
        ax : Matplotlib.Axes, optional
            Axes to plot the obs_points in, by default None which creates an Axes

        Returns
        -------
        axes: matplotlib axes instance

        """
        dfs = self.ml.read_obs_node()
        if ax is None:
            _, ax = plt.subplots(figsize=figsize, **kwargs)
        for i, df in dfs.items():
            df.plot(y=data, ax=ax, label=f"Node {i}", use_index=True)

        ax.set_xlabel(f"Time [{self.ml.basic_info['TUnit']}]")
        ax.set_ylabel(data)
        return ax
