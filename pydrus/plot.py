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
    
    def profile_information(self, data="Pressure Head", times = None, 
                             legend = True, figsize=(7, 5), 
                             title="Profile Information", 
                             nodes = "all", cmap="YlOrBr", **kwargs):
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
        title: str, optional
        cmap: str, optional
            String with a named Matplotlib colormap.

        Returns
        -------
        ax: matplotlib axes instance

        """
        l_unit = self.ml.basic_information["LUnit"]
        t_unit = self.ml.basic_information["TUnit"]  
        dfini = self.ml.read_nod_inf(times = 0)
      
        use_cols = ("Head","Moisture","K", "C", "Flux", "Sink")
        col_names = ("Pressure Head", "Water Content", 
                     "Hydraulic Conductivity", "Hydraulic Capacity", 
                     "Water Flux", "Root Uptake")
        units = ["h [{}]".format(l_unit),
                  "Theta [-]","K [{}/days]".format(l_unit),
                  "C [1/{}]".format(l_unit), 
                  "v [{}/{}]".format(l_unit, t_unit),
                  "S [1/{}]".format(t_unit)]
        
        col = col_names.index(data)        
        fig, ax = plt.subplots(figsize=figsize, **kwargs)
               
        if times is None:
            df = self.ml.read_nod_inf()
            for key, dataframe in df.items():
                ax.plot(dataframe[use_cols[col]], dataframe["Depth"],
                        label = "T" + str(key))
        else:  
            for time in times:
                df = self.ml.read_nod_inf(times = time)
                ax.plot(df[use_cols[col]], df["Depth"])
        
        space = (abs(dfini[use_cols[col]].min())-
                    abs(dfini[use_cols[col]].max()))        
        
        if data == "Pressure Head":
            ax.set_xlim(dfini[use_cols[col]].min(),
                        dfini[use_cols[col]].max()+space*0.2)
        if data == "Water Content":
            ax.set_xlim(dfini[use_cols[col]].min()+space*0.2,
                        dfini[use_cols[col]].max()-space*0.1)                      
        if data == "Hydraulic Conductivity":
            ax.set_xlim(0, dfini[use_cols[col]].max())            
        if data == "Hydraulic Capacity":
            ax.set_xlim(0)
        
        ax.plot(dfini[use_cols[col]], dfini["Depth"],
                    label = "T0", color = 'k')

        ax.set_ylim(self.ml.profile.loc[:, "x"].min(),
                    self.ml.profile.loc[:, "x"].max())
        ax.set_xlabel(units[col])
        ax.set_ylabel("Depth [{}]".format(self.ml.basic_information["LUnit"]))
        ax.set_title("Profile Information: " + data)
        ax.grid(linestyle='--')
        
        if legend:
            ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")    
        
        plt.tight_layout() 
        return ax
    
    def mass_balance(self, figsize=(6, 10), title="Mass_Balance_Information",
                        **kwargs):
        """Method to show the Mass balance information.

        Parameters
        ----------
        figsize: tuple, optional
        title: str, optional

        Returns
        -------
        balance: opens BALANCE.OUT in notepad
        """
        import subprocess
        
        folder = self.ml.ws_name
        balance = subprocess.Popen(["notepad.exe",folder + "\\BALANCE.OUT"])        
        return balance
    
    def water_flow(self, data="Potential Surface Flux", figsize=(10, 3), 
                   title="Water Flow", cmap="YlOrBr", **kwargs):
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
        title: str, optional
        cmap: str, optional
            String with a named Matplotlib colormap.

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
        
        cols = ("rTop", "rRoot", "vTop", "vRoot", "vBot", "hTop","hRoot", 
                "hBot", "RunOff", "Volume")
        df = self.ml.read_tlevel()
        col = col_names.index(data)
        
        if col < 5:
            fig, ax = plt.subplots(figsize=figsize, nrows = 1, ncols =2
                                   , **kwargs)
            ax[0].plot(df.index, df[cols[col]])
            ax[0].set_title(data)   
            ax[0].grid()            
            #Cumulative sum
            ax[1].plot(df.index, df["sum("+cols[col]+")"])
            ax[1].set_title("sum("+data+")")
            ax[1].grid()
        else:
            fig, ax = plt.subplots(figsize=(5,3), nrows = 1, ncols =1, 
                                   **kwargs)
            ax.plot(df.index, df[cols[col]])
            ax.set_title(data)   
            ax.grid() 
        return ax

    def shp(self, data="Water Content", figsize=(10, 3), 
            title="Soil hydraulic properties", cmap="YlOrBr", **kwargs):
        """Method to plot the soil hydraulic properties.

        Parameters
        ----------
        data: str, optional
            String with the variable of the water flow information to plot. 
            You can choose between: "Water Content", "Pressure head",
                     "log Pressure head", "Hydraulic Capacity",
                     "Hydraulic Conductivity", "log Hydraulic Conductivity",
                     "Effective Water Content".
            Default is "Water Content".        
        figsize: tuple, optional
        title: str, optional
        cmap: str, optional
            String with a named Matplotlib colormap.

        Returns
        -------
        ax: matplotlib axes instance

        """
        col_names = ("Water Content", "Pressure head",
                     "log Pressure head", "Hydraulic Capacity",
                     "Hydraulic Conductivity", "log Hydraulic Conductivity",
                     "Effective Water Content")
        cols = ("theta", "h", "log_h", "C", "K", "log_K", "S", "Kv")
        df = self.ml.read_I_check()
        col = col_names.index(data)

        fig, ax = plt.subplots(figsize=figsize, nrows = 1, ncols =3, **kwargs)
        fig.suptitle(data, fontsize=16, y=0.99)
        
        ax[0].plot(abs(df["h"]), df[cols[col]])  
        ax[0].grid()  
        ax[0].set_xlabel("h")         

        ax[1].plot(df["log_h"], df[cols[col]])  
        ax[1].grid()
        ax[1].set_xlabel("log_h")

        ax[2].plot(df["theta"], df[cols[col]])  
        ax[2].grid()
        ax[2].set_xlabel(r"$\dot{\Theta}$")
        return ax