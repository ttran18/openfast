#**********************************************************************************************************************************
# LICENSING
# Copyright (C) 2025 National Renewable Energy Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#**********************************************************************************************************************************

from ctypes import (
    CDLL,
    POINTER,
    create_string_buffer,
    byref,
    c_byte,
    c_int,
    c_double,
    c_float, 
    c_char,
    c_char_p, 
    c_bool,
    c_void_p
)
import numpy as np
import numpy.typing as npt
from pathlib import Path
import datetime
import os

from .interface_abc import OpenFASTInterfaceType

class SeaStateLib(OpenFASTInterfaceType):
    """
    This is the Python interface to the OpenFAST SeaState module.

    Notes:
    - SeaState is different from the other OpenFAST modules in that it does not do the typical
        CalcOutput-UpdateStates sequence. The calc_output function here is essentially an
        interrogation of the lookup table stored in SeaState. Therefore, the output_values
        attribute does not store a timeseries of values. It only stores a 1D array of outputs
        from the last call to calc_output.
    """

    def __init__(self, library_path):
        super().__init__(library_path)
        self.library_path = library_path

        self._initialize_routines()
        self.ended = False                  # For error handling at end

        # Create buffers for class data
        self.abort_error_level = 4
        self.error_status_c = c_int(0)
        self.error_message_c = create_string_buffer(self.ERROR_MSG_C_LEN)


        # This buffer for the channel names and units is set arbitrarily large
        # to start.  Channel name and unit lengths are currently hard
        # coded to 20 (this must match ChanLen in NWTC_Base.f90).
        self._channel_names_c = create_string_buffer(20 * 4000 + 1)
        self._channel_units_c = create_string_buffer(20 * 4000 + 1)

        self.numResPts = 0                  # Number of wind points we will
                                            # request information from
                                            # non-CalcOutput routines.

        self.numChannels = 0                # Number of channels returned

        # flags
        self.debuglevel  = 0                # 0-4 levels


    def _initialize_routines(self):
        self.SeaSt_C_PreInit.argtypes = [
            POINTER(c_float),       # intent(in   ) :: Gravity_c
            POINTER(c_float),       # intent(in   ) :: WtrDens_c
            POINTER(c_float),       # intent(in   ) :: WtrDpth_c
            POINTER(c_float),       # intent(in   ) :: MSL2SWL_c
            POINTER(c_int),         # intent(in   ) :: debuglevel
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char),        # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.SeaSt_C_PreInit.restype = c_int

        self.SeaSt_C_Init.argtypes = [
            POINTER(c_char_p),      # intent(in   ) :: InputFile_c(IntfStrLen)
            POINTER(c_char_p),      # intent(in   ) :: OutRootName_c(IntfStrLen)
            POINTER(c_int),         # intent(in   ) :: NSteps_c
            POINTER(c_float),       # intent(in   ) :: TimeInterval_c
            POINTER(c_int),         # intent(in   ) :: WaveElevSeriesFlag_c
            POINTER(c_int),         # intent(in   ) :: WrWvKinMod_c
            POINTER(c_int),         # intent(  out) :: NumChannels_c
            POINTER(c_char),        # intent(  out) :: OutputChannelNames_C
            POINTER(c_char),        # intent(  out) :: OutputChannelUnits_C
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char),        # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.SeaSt_C_Init.restype = c_int

        self.SeaSt_C_CalcOutput.argtypes = [
            POINTER(c_double),      # intent(in   ) :: Time_C
            POINTER(c_float),       # intent(  out) :: OutputChannelValues_C(p%NumOuts)
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char),        # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.SeaSt_C_CalcOutput.restype = c_int

        self.SeaSt_C_End.argtypes = [
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char)         # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]
        self.SeaSt_C_End.restype = c_int

        self.SeaSt_C_GetWaveFieldPointer.argtypes = [
            POINTER(c_void_p),      # intent(  out) :: pointer to the WaveField data
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char),        # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]

        self.SeaSt_C_SetWaveFieldPointer.argtypes = [
            POINTER(c_void_p),      # intent(in   ) :: pointer to the WaveField data
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char),        # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]


        self.SeaSt_C_GetFluidVelAccDens.argtypes = [
            POINTER(c_double),      # intent(in   ) :: Time_C
            POINTER(c_float),       # intent(in   ) :: Pos_c(3)
            POINTER(c_float),       # intent(  out) :: Vel_c(3)
            POINTER(c_float),       # intent(  out) :: Acc_c(3)
            POINTER(c_int),         # intent(  out) :: NodeInWater_C
            POINTER(c_float),       # intent(  out) :: Density_C
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char)         # intent(  out) :: ErrMsg_C(ErrMsgLen_C)

        ]

        self.SeaSt_C_GetSurfElev.argtypes = [
            POINTER(c_double),      # intent(in   ) :: Time_C
            POINTER(c_float),       # intent(in   ) :: Pos_c(3)
            POINTER(c_float),       # intent(  out) :: Elev_C
            POINTER(c_int),         # intent(  out) :: ErrStat_C
            POINTER(c_char)         # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
        ]

    def check_error(self) -> None:
        """Checks for and handles any errors from the Fortran library.

        Raises:
            RuntimeError: If a fatal error occurs in the Fortran code
        """
        # If the error status is 0, return
        if self.error_status_c.value == 0:
            return

        # Get the error level and error message
        error_level = self.error_levels.get(
            self.error_status_c.value,
            f"Unknown Error Level: {self.error_status_c.value}"
        )
        error_msg = self.error_message_c.raw.decode('utf-8').strip()
        message = f"WaveTank library {error_level}: {error_msg}"
        # If the error level is fatal, call WaveTank_End() and raise an error
        if self.error_status_c.value >= self.abort_error_level:
            try:
                self.SeaSt_C_End(
                    byref(self.error_status_c),             # OUT <- error status code
                    self.error_message_c                    # OUT <- error message buffer
                )
                if self.error_status_c.value == 4:
                    error_msg = self.error_message_c.raw.decode('utf-8').strip()
                    print(f'WaveTank_End error: {error_msg}')
            except Exception as e:
                message += f"\nAdditional error during cleanup: {e}"
            raise RuntimeError(message)
        else:
            print(message)


    #FIXME: store these elsewhere
    def seastate_preinit(
        self,
        gravity: float = 9.80665,
        water_density: float = 1025,
        water_depth: float = 200,
        msl2swl: float = 0,
    ):
        """Set environment variables and general setup

        Args:

        Raises:
            ValueError: If values are outside reasonable bounds
            RuntimeError: If preinit fails 
        """
        self.SeaSt_C_PreInit(
            byref(c_float(gravity)),
            byref(c_float(water_density)),
            byref(c_float(water_depth)),
            byref(c_float(msl2swl)),
            byref(c_int(self.debug_level)),         # IN -> debug level (0=None to 4=all meshes)
            byref(self.error_status_c),             # OUT <- error status code
            self.error_message_c                    # OUT <- error message buffer
        )
        self.check_error()


    def seastate_init(
        self,
        primary_ss_file,
        outrootname: str = "./seastate.SeaSt",
        wave_kinematics_mode: int = 0,
        n_steps: int = 801,
        time_interval: float = 0.125,
        wave_elevation_series_flag: int = 0,
    ):

        # This buffer for the channel names and units is set arbitrarily large
        # to start. Channel name and unit lengths are currently hard
        # coded to 20 (this must match ChanLen in NWTC_Base.f90).
        _channel_names = create_string_buffer(20 * 4000 + 1)
        _channel_units = create_string_buffer(20 * 4000 + 1)
        self._numChannels = c_int(0)


        self.SeaSt_C_Init(
            c_char_p(primary_ss_file.encode('utf-8')),
            c_char_p(outrootname.encode('utf-8')),
            byref(c_int(n_steps)),
            byref(c_float(time_interval)),
            byref(c_int(wave_elevation_series_flag)),
            byref(c_int(wave_kinematics_mode)),
            byref(self._numChannels),
            _channel_names,
            _channel_units,
            byref(self.error_status_c),             # OUT <- error status code
            self.error_message_c                    # OUT <- error message buffer
        )
        self.check_error()

        # Initialize output channels
        self.numChannels = self._numChannels.value

        # if len(_channel_names.value.split()) == 0:
        #     self.output_channel_names = []
        # else:
        #     self.output_channel_names = [n.decode('UTF-8') for n in _channel_names.value.split()] 
        self.output_channel_names = [n.decode('UTF-8') for n in _channel_names.value.split()] 

        # if len(_channel_units.value.split()) == 0:
        #     self.output_channel_units = []
        # else:
        #     self.output_channel_units = [n.decode('UTF-8') for n in _channel_units.value.split()] 
        self.output_channel_units = [n.decode('UTF-8') for n in _channel_units.value.split()] 

        # Allocate the data for the outputs
        self.output_values = np.zeros( self._numChannels.value, dtype=c_float, order='C' )

    def seastate_calcOutput(self, time: float, output_channel_values: npt.NDArray[np.float32]) -> None:
        """Calculate output values at the given time.

        Args:
            time: Current simulation time
            output_channel_values: Array to store calculated output values

        Raises:
            ValueError: If output_channel_values array has wrong size
            RuntimeError: If calculation fails
        """
        if output_channel_values.size != self.numChannels:
            raise ValueError(
                f"Output array must have size {self.numChannels}, "
                f"got {output_channel_values.size}"
            )

        output_channel_values_c = (c_float * self.numChannels)(0.)

        self.SeaSt_C_CalcOutput(
            byref(c_double(time)),           # IN -> current simulation time
            self.output_values.ctypes.data_as(POINTER(c_float)), # OUT: output channel values
            byref(self.error_status_c),      # OUT <- error status
            self.error_message_c             # OUT <- error message
        )
        self.check_error()

        # Copy results back to numpy array
        output_channel_values[:] = np.reshape(self.output_values, (self.numChannels))

    def seastate_end(self):
        if not self.ended:
            self.ended = True

            self.SeaSt_C_End(
                byref(self.error_status_c),             # OUT <- error status code
                self.error_message_c                    # OUT <- error message buffer
            )

        self.check_error()



    def seastate_getWaveFieldPointer(self,ss_pointer: c_void_p) -> None:
        self.SeaSt_C_GetWaveFieldPointer(
            byref(ss_pointer),                          # IN  -> pointer to the WaveField data
            byref(self.error_status_c),                 # OUT <- error status code
            self.error_message_c                        # OUT <- error message buffer
        )
        self.check_error()

    def seastate_setWaveFieldPointer(self,ss_pointer: c_void_p) -> None:
        self.SeaSt_C_SetWaveFieldPointer(
            byref(ss_pointer),                          # IN  -> pointer to the WaveField data
            byref(self.error_status_c),                 # OUT <- error status code
            self.error_message_c                        # OUT <- error message buffer
        )
        self.check_error()


#    def get_fluidVelAccDens(self, pos, t):
#
#SeaSt_C_GetFluidVelAccDens.argtypes = [
#            POINTER(c_double),      # intent(in   ) :: Time_C
#            POINTER(c_float),       # intent(in   ) :: Pos_c(3)
#            POINTER(c_float),       # intent(  out) :: Vel_c(3)
#            POINTER(c_float),       # intent(  out) :: Acc_c(3)
#            POINTER(c_int),         # intent(  out) :: NodeInWater_C
#            POINTER(c_float),       # intent(  out) :: Density_C
#            POINTER(c_int),         # intent(  out) :: ErrStat_C
#            POINTER(c_char)         # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
#
#        ]
#
#        self.SeaSt_C_GetSurfElev.argtypes = [
#            POINTER(c_double),      # intent(in   ) :: Time_C
#            POINTER(c_float),       # intent(in   ) :: Pos_c(3)
#            POINTER(c_float),       # intent(  out) :: Elev_C
#            POINTER(c_int),         # intent(  out) :: ErrStat_C
#            POINTER(c_char)         # intent(  out) :: ErrMsg_C(ErrMsgLen_C)
#        ]


    @property
    def num_outs(self):
        return self._numChannels.value


#===============================================================================
#   Helper classes for writing output channels to file.
#   For the regression testing to mirror the output from the InfowWind Fortran
#   driver.  This may also have value for debugging the interfacing to SS.

class ResultsOut():
    """
    This is only for testing purposes. Since we are not returning the
    velocities to anything, we will write them to file as we go for
    comparison in the regression test.  When coupled to another code, the
    velocities array would be passed back to the calling code for use in
    the aerodynamic solver.
    """
    def __init__(self, filename, NumResPts):

        self.results_file = open(filename, 'w')        # open output file and write header info

        # write file header
        t_string=datetime.datetime.now()
        dt_string=datetime.date.today()
        self.results_file.write(f"## This file was generated by SeaState called from Python on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}{os.linesep}")
        self.results_file.write(f"## This file contains outputs from calls to SeaState routines (not the CalcOutput) at the {NumResPts} points specified in the file {filename}{os.linesep}")
        self.results_file.write(f"# {os.linesep}")
        self.results_file.write(f"# {os.linesep}")
        self.results_file.write(f"# {os.linesep}")
        self.results_file.write(f"# {os.linesep}")
        self.results_file.write(f"          T                  X                  Y                  Z                  U                  V                  W{os.linesep}")
        self.results_file.write(f"         (s)                (m)                (m)                (m)               (m/s)              (m/s)              (m/s){os.linesep}")
        self.opened = True

    def write(self,t,positions,velocities):
        for p, v in zip(positions,velocities):
            # TODO: does \n work as expected on Windows?
            self.results_file.write('  %16.8f   %16.8f   %16.8f   %16.8f   %16.8f   %16.8f   %16.8f\n' % (t,p[0],p[1],p[2],v[0],v[1],v[2]))

    def end(self):
        if self.opened:
            self.results_file.close()
            self.opened = False




class WriteOutChans():
    """
    This is only for testing purposes. Since we are not returning the
    output channels to anything, we will write them to file.  When coupled to
    another code, this data would be passed back for inclusion the any output
    file there.
    """
    def __init__(self,filename,chan_names,chan_units):
        chan_names.insert(0,'Time')             # add time index header
        chan_units.insert(0,'(s)')              # add time index unit
        self.OutFile=open(filename,'wt')        # open output file and write header info
        # write file header
        t_string=datetime.datetime.now()
        dt_string=datetime.date.today()
        self.OutFile.write(f"## This file was generated by SeaState c-bindings library on {dt_string.strftime('%b-%d-%Y')} at {t_string.strftime('%H:%M:%S')}\n")
        self.OutFile.write(f"## This file contains output channels requested from the OutList section of the input file")
        self.OutFile.write(f"{filename}\n")
        self.OutFile.write("#\n")
        self.OutFile.write("#\n")
        self.OutFile.write("#\n")
        self.OutFile.write("#\n")
        l = len(chan_names)
        f_string = "{:^15s}"+"   {:^20s}  "*(l-1)
        self.OutFile.write(f_string.format(*chan_names) + '\n')
        self.OutFile.write(f_string.format(*chan_units) + '\n')
        self.opened = True

    def write(self,chan_data):
        l = chan_data.shape[1]
        f_string = "{:10.4f}"+"{:25.7f}"*(l-1)
        for i in range(0,chan_data.shape[0]):
            self.OutFile.write(f_string.format(*chan_data[i,:]) + '\n')
            #if i==0:
            #    print(f"{chan_data[i,:]}")

    def end(self):
        if self.opened:
            self.OutFile.close()
            self.opened = False
