
"""
xltek_converter.py

Last Edited: 12/1/2017

Lead Author[s]: Anthony Fong
Contributor[s]:

Description:


Machine I/O
Input:
Output:

User I/O
Input:
Output:


"""
########################################################################################################################

########## Libraries, Imports, & Setup ##########

# Default Libraries #
import os
import pathlib
import datetime
import multiprocessing

# Downloaded Libraries #
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Imports from Local Packages #
import xltek_reader




########## Definitions ##########

# Classes #

class XLTEK_Converter:
    def __init__(self, subj, subs_dir=pathlib.Path('Z:/human/converted_clinical'), reader=None):
        ## Attributes ##
        self.subj = subj
        if isinstance(subs_dir, pathlib.Path):
            self.subs_path = subs_dir
        else:
            self.subs_path = pathlib.Path(subs_dir)
        self.subj_dir = pathlib.Path(self.subs_path, self.subj)
        if reader is None:
            self.reader = xltek_reader.XLTEK_Reader(subj, update=False)
        else:
            self.reader = reader
        self.day_dirs = []
        self.h5_fobj = None
        self.pos = (0, 0, 0, -1)
        self.cday_index = 0
        self.current_day = None
        self.current_file = None

        self.stop_events = {'convert':multiprocessing.Event()}

        if not self.subj_dir.is_dir():
            print('No such subject directory as %s'%(str(self.subj_dir)))
            ans = input('Would you like to create a new subject directory? (y/N): ')
            if ans.upper() == 'Y' or ans.upper() == 'Yes':
                self.subj_dir.mkdir()

        self.build_directories()
        self.current_day = self.day_dirs[0]
        self.aquire_conv_point()

    def __getstate__(self):
        state = self.__dict__.copy()
        state['reader'] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def build_directories(self):
        start = self.reader._studies[0].creation_time.date()
        end = self.reader._studies[-1].end_time.date()
        days = end - start
        for i in range(0,days.days+1):
            day = start + datetime.timedelta(days=i)
            self.day_dirs.append(Day_Directory(self.subj, day, self.subj_dir))

    def next_day(self, new=None):
        index = self.cday_index + 1
        if index >= len(self.day_dirs):
            if isinstance(new, datetime.date):
                self.day_dirs.append(Day_Directory(self.subj, new, self.subj_dir))
            else:
                index = 0
        self.cday_index = index
        return self.day_dirs[index]

    def aquire_conv_point(self):
        for i in range(len(self.day_dirs)-1, -1, -1):
            if len(self.day_dirs[i].data_files) > 0:
                self.current_day = self.day_dirs[i]
                self.current_file = self.current_day.data_files[-1]
                self.current_file.open_h5()
                self.cday_index = i
                self.pos = self.current_file.end_entry
                break
        return self.pos

    def spawn_collector(self, start_entry):
        self.buffer, self.forward_process, self.stop_events['buffer'] = self.reader.start_forward(start_entry, spawn=True)
        return self.buffer

    def format_entry(self, entry):
        data = entry['data']
        dshape = data.shape
        channels = np.arange(0, dshape[1])
        # print('snc: %d start: %d'%(entry['snc_sample'],entry['start_sample']))
        samples = np.array(range(entry['start_sample'], entry['end_sample']))

        locs = np.zeros((dshape[0],4),dtype=np.int32)
        locs[:,:] = entry['entry_info']
        #print(entry['snc_sample'])

        times = np.zeros(dshape[0], dtype = np.float64)
        for i in range(0, len(samples)):
            delta_t = datetime.timedelta(seconds=((samples[i] - entry['snc_sample']) * 1.0 / entry['sample_rate']))
            time = entry['snc_time'] + delta_t
            times[i] = time.timestamp()

        return data, dshape, samples, times, locs, channels, entry['sample_rate']

    def parse_entry(self, entry):
        # print(entry['entry_info'])
        if entry['entry_info'][:3] == self.pos[:3]:
            return
        elif entry['entry_info'][:2] == self.pos[:2]:
            self.current_file.add_data(entry)
        else:
            # print('new file')
            self.current_file.close_h5()
            ent_time = entry['snc_start']
            if ent_time.date() > self.current_day.day:
                self.current_day = self.next_day(ent_time.date())
            self.current_file = self.current_day.new_file(entry)
        self.pos = entry['entry_info']

    def start_convert(self):
        self.spawn_collector(self.pos)
        self.convert_process = multiprocessing.Process(target=self._con_process, name='Convert Process')
        self.convert_process.daemon = True
        self.convert_process.start()

    def stop_convert(self):
        self.stop_events['convert'].set()
        self.stop_events['buffer'].set()
        del self.buffer[:]
        self.convert_process.join()
        self.forward_process.join()
        del self.buffer[:]
        self.stop_events['convert'].clear()
        self.stop_events['buffer'].clear()

    def _con_process(self, buffer=None):
        if buffer is None:
            buffer = self.buffer
        while not self.stop_events['convert'].is_set() and self.current_file is None:
            if len(buffer) > 0:
                entry = buffer.pop(0)
                self.current_file = self.current_day.new_file(entry)
                self.pos = entry['entry_info']
        while not self.stop_events['convert'].is_set():
            if len(buffer) > 0:
                # print(len(buffer))
                self.parse_entry(buffer.pop(0))


class Day_Directory:
    def __init__(self, subj, day, subj_dir, path=None):
        self.subj = subj
        self.day = day

        if isinstance(subj_dir, pathlib.Path):
            self.subj_dir = subj_dir
        else:
            self.subj_dir = pathlib.Path(subj_dir)

        if isinstance(path, pathlib.Path):
            self.path = path
        elif isinstance(path, str):
            self.path = pathlib.Path(path)
        else:
            self.path = pathlib.Path(self.subj_dir, subj+'_'+str(day))

        self.data_files = []

        if not self.path.is_dir():
            self.path.mkdir()
        else:
            self.data_files = [Data_File(self.subj, path=p) for p in self.path.glob('*.h5')]

    def new_file(self, entry):
        file = Data_File(self.subj, self.path, entry)
        self.data_files.append(file)
        return file

class Data_File:
    def __init__(self, subj, path=None, entry=None):
        self.subj = subj
        if isinstance(path, pathlib.Path):
            self.path = path
        else:
            self.path = pathlib.Path(path)
        self.name = path.parts[-1]
        self.start = None
        self.is_open = True
        self.start_entry = None
        self.end_entry = None

        if self.path.is_file():
            self.open_h5()
            self.close_h5()
        elif self.path.is_dir() and entry is not None:
            self.dir = self.path
            self.start = entry['snc_start'].replace(tzinfo=None)
            self.start_entry = entry['entry_info']
            self.end_entry = entry['entry_info']
            self.name = subj + '_' + self.start.isoformat('_','seconds').replace(':','~')
            self.path = pathlib.Path(path, self.name+'.h5')
            if self.path.is_file():
                self.name = subj + '_' + self.start.isoformat('_', 'milliseconds').replace(':', '~')
                self.path = pathlib.Path(path, self.name + '.h5')
                if self.path.is_file():
                    self.name = subj + '_' + self.start.isoformat('_', 'microseconds').replace(':', '~')
                    self.path = pathlib.Path(path, self.name + '.h5')
            self.make_h5file(self.path, entry)

    def __repr__(self):
        return repr((self.start))

    def __getstate__(self):
        state = self.__dict__.copy()
        name = str(state['path'])
        open = state['is_open']
        if open:
            fobj = state['h5_fobj']
            fobj.flush()
            fobj.close()
        state['h5_fobj'] = (name, open)
        return state

    def __setstate__(self, state):
        name, open = state['h5_fobj']
        state['h5_fobj'] = h5py.File(str(name),'r+')
        if not open:
            state['h5_fobj'].close()
        self.__dict__.update(state)

    def format_entry(self, entry):
        data = entry['data']
        dshape = data.shape
        channels = np.arange(0, dshape[1])
        # print('snc: %d start: %d'%(entry['snc_sample'],entry['start_sample']))
        samples = np.array(range(entry['start_sample'], entry['end_sample']))

        locs = np.zeros((dshape[0],4),dtype=np.int32)
        locs[:,:] = entry['entry_info']
        #print(entry['snc_sample'])

        times = np.zeros(dshape[0], dtype = np.float64)
        for i in range(0, len(samples)):
            delta_t = datetime.timedelta(seconds=((samples[i] - entry['snc_sample']) * 1.0 / entry['sample_rate']))
            time = entry['snc_time'] + delta_t
            times[i] = time.timestamp()

        return data, dshape, samples, times, locs, channels, entry['sample_rate']

    def make_h5file(self, f_path, entry):
        data, dshape, samples, times, locs, channels, sample_rate = self.format_entry(entry)

        self.h5_fobj = h5py.File(str(f_path))
        self.h5_fobj.attrs['start time'] = times[0]
        self.h5_fobj.attrs['end time'] = times[-1]
        self.h5_fobj.attrs['total samples'] = len(samples)
        self.h5_fobj.attrs['start entry'] = locs[0]
        self.h5_fobj.attrs['end entry'] = locs[0]

        ecog = self.h5_fobj.create_dataset('ECoG Array', dtype='f8', data=data, maxshape=(None,None), compression="gzip", compression_opts=4)
        ecog.attrs['Sampling Rate'] = sample_rate

        samplestamps = self.h5_fobj.create_dataset('samplestamp vector', dtype='i', data=samples, maxshape=(None,), compression="gzip", compression_opts=4)
        timestamps = self.h5_fobj.create_dataset('timestamp vector', dtype='f8', data=times, maxshape=(None,), compression="gzip", compression_opts=4)
        entrystamps = self.h5_fobj.create_dataset('entry vector', dtype='i', data=locs, maxshape=(None,4), compression="gzip", compression_opts=4)
        channelstamps = self.h5_fobj.create_dataset('channel indices', dtype='i', data=channels, maxshape=(None,), compression="gzip", compression_opts=4)

        ecog.dims.create_scale(samplestamps, 'sample axis')
        ecog.dims.create_scale(timestamps, 'time axis')
        ecog.dims.create_scale(entrystamps, 'entry axis')
        ecog.dims.create_scale(channelstamps, 'channel axis')

        ecog.dims[0].attach_scale(samplestamps)
        ecog.dims[0].attach_scale(timestamps)
        ecog.dims[0].attach_scale(entrystamps)
        ecog.dims[1].attach_scale(channelstamps)

        self.h5_fobj.flush()

        return self.h5_fobj

    def add_data(self, entry, h5_fobj=None):
        if h5_fobj is None:
            h5_fobj = self.h5_fobj

        data, dshape, samples, times, locs, channels, sample_rate = self.format_entry(entry)

        self.end_entry = tuple(locs[0])
        h5_fobj.attrs['end entry'] = locs[0]
        h5_fobj.attrs['end time'] = times[-1]

        ecog = h5_fobj['ECoG Array']
        last_total = ecog.shape[0]

        samplestamps = h5_fobj['samplestamp vector']
        timestamps = h5_fobj['timestamp vector']
        entrystamps = h5_fobj['entry vector']

        ecog.resize(last_total+dshape[0], 0)
        samplestamps.resize(last_total+dshape[0], 0)
        timestamps.resize(last_total+dshape[0], 0)
        entrystamps.resize(last_total+dshape[0], 0)

        ecog[-dshape[0]:,:] = data
        samplestamps[-dshape[0]:] = samples
        timestamps[-dshape[0]:] = times
        entrystamps[-dshape[0]:,:] = locs

        h5_fobj.attrs['total samples'] = len(samplestamps)

    def open_h5(self):
        try:
            self.h5_fobj = h5py.File(str(self.path))
            self.start = datetime.datetime.fromtimestamp(self.h5_fobj.attrs['start time'])
            self.start_entry = tuple(self.h5_fobj.attrs['start entry'])
            self.end_entry = tuple(self.h5_fobj.attrs['end entry'])
            self.is_open = True
        except:
            print(str(self.path))
            #print(self.h5_fobj.attrs['start time'])
            print(datetime.datetime.fromtimestamp(self.h5_fobj.attrs['start time']))
            raise NameError('Something went wrong')

    def close_h5(self):
        self.h5_fobj.flush()
        self.h5_fobj.close()
        self.is_open = False


if __name__ == "__main__":
    bland = XLTEK_Converter('EC164')