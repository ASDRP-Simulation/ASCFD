import configparser
import numpy as np
import os


class CaseInsensitiveConfigParser(configparser.ConfigParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.optionxform = (
            str.lower
        )  # Transform options to lowercase for case insensitivity


class Inputs:
    def __init__(self, fname):
        config = CaseInsensitiveConfigParser()
        config.read(fname)

        if not os.path.exists(fname):
            raise RuntimeError("Inputs File not found.")

        # Mesh
        self.nx = self.get_config_value(config, "Mesh", "nx", type_func=int)
        self.numghosts = self.get_config_value(
            config, "Mesh", "numghosts", type_func=int
        )
        self.x_lo = self.get_config_value(config, "Mesh", "x_lo", type_func=float)
        self.x_hi = self.get_config_value(config, "Mesh", "x_hi", type_func=float)

        self.xlim = (self.x_lo, self.x_hi)

        # Time
        self.nt = self.get_config_value(
            config, "Time", "time_steps", type_func=int, mandatory=False, default=np.inf
        )
        self.timeStepper = self.get_config_value(config, "Time", "method")
        self.t0 = self.get_config_value(
            config, "Time", "t0", mandatory=False, default=0.0, type_func=float
        )
        self.t_finish = self.get_config_value(
            config, "Time", "t_finish", mandatory=False, default=np.inf, type_func=float
        )

        if (self.t_finish == np.inf) and (self.time_steps == np.inf):
            raise RuntimeError(
                "Must supply either t_finish or time_steps in [Time] section of inputs file."
            )
        self.cfl = self.get_config_value(config, "Time", "cfl", type_func=float)


        if self.t0 >= self.t_finish:
            raise RuntimeError("Initial time is >= to the final time.")

        # Fluid
        self.ics = self.get_config_value(config, "Fluid", "ics")
        self.system = self.get_config_value(config, "Fluid", "system")

        #if self.system == "euler1DNS2":
        self.gammas = self.get_config_value(config, "Fluid", "gammas", type_func=lambda s: [float(item) for item in s.strip('[]').split(',')])
        self.mW = self.get_config_value(config, "Fluid", "mW", type_func=lambda s: [float(item) for item in s.strip('[]').split(',')])

        # Method
        self.flux = self.get_config_value(config, "Method", "flux")
        self.bc_lo = self.get_config_value(config, "Method", "bcs_lo")
        self.bc_hi = self.get_config_value(config, "Method", "bcs_hi")

        # Output
        self.output_freq = self.get_config_value(
            config, "Output", "output_freq", type_func=int, mandatory=False, default=1
        )
        self.output_dir = self.get_config_value(
            config, "Output", "output_dir", mandatory=False, default="output/"
        )
        self.make_movie = self.get_config_value(
            config,
            "Output",
            "make_movie",
            type_func=bool,
            mandatory=False,
            default=True,
        )


    def get_config_value(
        self, config, section, option, type_func=str, mandatory=True, default=None
    ):
        try:
            # Attempt to get the value with type conversion; if not found or conversion fails, it will raise an error
            value = config.get(
                section, option, fallback=default if not mandatory else None
            )
            if value is None and mandatory:
                raise ValueError(
                    f"Missing mandatory argument: '{option}' in section '{section}'"
                )
            return type_func(value)
        except (configparser.NoSectionError, configparser.NoOptionError) as e:
            if mandatory:
                raise ValueError(
                    f"Missing mandatory argument: '{option}' in section '{section}'"
                ) from e
            return default
        except ValueError as e:
            raise ValueError(
                f"Type conversion error for '{option}' in section '{section}': {e}"
            ) from e
