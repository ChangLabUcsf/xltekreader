
"""
xltek_reader.py

Last Edited: 12/1/2017

Lead Author[s]: Anthony Fong, Pier Mantovani
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
import sys, os, pathlib
import struct
import re

import multiprocessing
import multiprocessing.pool
import time, datetime, pytz

# Downloaded Libraries #
# import pymysql
# pymysql.install_as_MySQLdb()
import numpy as np

# Imports from Local Packages #
# cur_dir = os.getcwd()
# sys.path.append(os.path.join(cur_dir, 'cldb-framework/changlabdb'))
# from changlabdb.backend import *
from xlrdlib2.x64.win32_1_0_rc1 import libpyxlrd



########## Definitions ##########

_FILETIME_null_date = datetime.datetime(1601, 1, 1, 0, 0, 0, tzinfo=pytz.utc)
ucsf_time = pytz.timezone('America/Los_Angeles')


# Classes #
class XLTEK_Reader:
    """
    XLTEK_reader: An object that contains contains all patient's studies created by the XLTEK system.

    Required Modules: os, pathlib, struct, re, multiprocessing, time, datetime, pytz, changlabdb, xlrdlib2, numpy
    Required Classes: XLTEK_Study
    Methods: get_XLTEK_ID, get_neuro_dirs, get_subj_name, get_studies, create_studies, get_total_entries, last_ent,
             find_entry, get_data_range, get_last_data, load_buffer, update_study_dirs, update_studies, start_updater,
             stop_updater, start_current, stop_current, start_forward, stop_all_forward, shutdown, _studies_updater,
             _buffer_process, _forward_process, _current_process, proxy, proxies

    Class Attributes
    none

    Object Parameters & Attributes
    Parameters
    :param subj:        string      The subject's common ID name.  
    :param study_ID:    string      The XLTEK ID corresponding to the first study.  
    :param network_dir: string      The directory that contains all the neuroworks network folders. 
    :param update:      boolean     A flag to start updating the buffer after initializing the XLTEK_Reader object.
    :param manager:     object      Optional multiprocessing.Mananger object that will hold all of the XLTEK_Reader attributes.
    
    Attributes
    :subj:              string      The subject's common ID name.
    :f_study_ID:        string      The XLTEK ID corresponding to the first study. 
    :_first_name:       string      The first name of the subject.
    :_last_name:        string      The last name of the subject.
    :_name_key:         string      The name key used by XLTEK to identify a subjects studies.
    
    :beginning:         datetime    The beginning of the first study.
    :end_time:          datetime    The end of the last study.
    :buffer_limit:      double      The number of data entries that the buffer will contain.
    :current_study:     double      The index of the most recent study.
    :total_data_entries:double      The total number of data entries for this subject. 
    
    :network_dir:       path        The directory that contains all the neuroworks network folders.
    :neuroworks_dirs:   list        A Manager list of paths the neuroworks directory paths.
    :_all_study_dirs:   set         A set of paths of all the studies for the subject even the irrelevant ones. 
    :studies:           list        A Manager list of relevant XLTEK_Studies of the subject.
    :_studies:          list        A normal list of the XLTEK_Studies meant for un-updated quick accesses.
    :_buffer:           list        A Manager list that will contain streaming data.
    :_for_buffer:       list        A Manager list that will contain data requested from the forward process.

    :manager:           object      The multiprocessing.Mananger object that maintains all of the XLTEK_Reader attributes.
    :proxy_dict:        dictionary  A Manager dictionary that holds all of the XLTEK_Reader attributes.
    :stop_events:       dictionary  A dictionary containing all of the stop Events objects for multiprocesses.

    :studies_process:   object      A process that checks if there are new studies and adds them to the studies list
    :current_process:   object      A process that keeps that data up to date with the past hour's worth of data
    :forward_process:   list        A list of processes created to return data from a given time to the most recent data
    """
    def __getstate__(self):
        state = self.__dict__.copy()
        state['manager'] = None
        return state

    def __setstate__(self, state):
        names = dict(state['proxy_dict'])
        for key in names:
            self.proxy(key)
        self.__dict__.update(state)

    def __del__(self):
        self.shutdown()

    def __init__(self, subj, study_ID=None, neuro_dirs=None, network_dir=pathlib.Path('//mb-files.ucsfmedicalcenter.org/neuroworksshare$'), update=True, manager=None):
        ## Proxy Manager Setup ##
        if manager and isinstance(manager, multiprocessing.Manager):
            self.manager = manager
        else:
            self.manager = multiprocessing.Manager()
        self.proxy_dict = self.manager.dict()

        ## Proxy Attributes ##
        # Identifiers #
        self.subj = subj
        if study_ID is None:
            # self.f_study_ID = self.get_XLTEK_ID(subj)
            raise ValueError("Missing Subject ID")
        else:
            self.f_study_ID = study_ID
        self._first_name = ''
        self._last_name = ''
        self._name_key = ''

        # Properties #
        self.beginning = None

        self.buffer_limit = 3600
        self.current_study = 0

        self.network_dir = pathlib.Path(network_dir)
        self._all_study_dirs = set()

        self.proxies(exclude=('proxy_dict','manager'))

        # Proxy Data Arrays #
        self.neuroworks_dirs = self.manager.list()
        self.studies = self.manager.list()
        self._buffer = self.manager.list()
        self._for_buffer = self.manager.list()

        ## Local Attributes ##
        self._end_time = None

        # Data Arrays #
        self._studies = []

        # Multiprocessing  #
        self.stop_events = {'update': multiprocessing.Event(), 'fill': multiprocessing.Event(),
                            'forward':[], 'current':multiprocessing.Event()}

        ## Initialization Methods ##
        self.get_neuro_dirs(neuro_dirs)
        self.get_subj_name(self.f_study_ID)
        self.get_studies()
        self._studies = list(self.studies)

        # Initialize Multiprocessing Objects #
        self.studies_process = multiprocessing.Process(target=self._studies_updater, name='Update Studies')
        self.current_process = multiprocessing.Process(target=self._current_process, name='Current Process',
                                                       args=(self.stop_events['current'], self._buffer, True))
        self.forward_process = []

        if update:
            self.start_updater()


    # def get_XLTEK_ID(self, subj):
    #     '''
    #     get_XLTEK_ID: Uses the database to return the XLTEK ID
    #
    #     Parameters
    #     :param subj:    string      The subject lab reference ID
    #
    #     Returns
    #     :xltek_num:        string      The XLTEK ID
    #     '''
    #     xltek_num = 0
    #     for i in PatientPropertyRecord.property_type.get_queryset(PatientPropertyDefinition="First Xltek Study ID"):
    #         if "First Xltek Study ID" in str(i):
    #             for j in i.patientpropertyrecord_set.all():
    #                 if str(subj) in str(j.patient):
    #                     xltek_num = j.property_value
    #                     return xltek_num
    #     if xltek_num == 0:
    #         raise NameError("Subject %s does not have XLTEK Number inputted in Database, please input the XLTEK number "
    #                         "in database and try again." %(subj))

    def get_neuro_dirs(self, neuroworks_dirs=None):
        '''
        get_neruo_dirs: Assigns neuroworks_dirs as a list of all Neuroworks directories in the Network folder. 
        '''
        print('Finding Network Directories...')
        if neuroworks_dirs is None:
            for location in self.network_dir.glob('Neuroworks*'):
                try:
                    sub_dir = pathlib.Path(location, 'DbData')
                    if sub_dir.is_dir():
                        os.path.getatime(str(sub_dir))
                        self.neuroworks_dirs.append(sub_dir)
                except:
                    print('Cannot access %s skipping.' % (location))
        else:
            new_dirs = []
            if isinstance(neuroworks_dirs, pathlib.Path):
                new_dirs = [neuroworks_dirs]
            elif neuroworks_dirs is not list:
                neuroworks_dirs = [neuroworks_dirs]
            for dir_ in neuroworks_dirs:
                if not isinstance(dir_, pathlib.Path):
                    dir_ = pathlib.Path(dir_)
                new_dirs.append(dir_)
            self.neuroworks_dirs = new_dirs

        print('%d Network Directories Found.'%(len(self.neuroworks_dirs)))

    def get_subj_name(self, study_ID):
        '''
        get_subj_name: Gets the subject's name and information about the first study, assigning this XLTEK_Reader's attributes.
        
        Parameters
        :param study_ID:    string  The XLTEK study ID
        '''
        print('Finding Subject %s...'%(self.subj))
        for location in self.neuroworks_dirs:
            study_dirs = location.glob(('*'+study_ID))
            for study_dir in study_dirs:
                # Get Name
                study_name = str(os.path.split(study_dir)[-1])
                last_end = study_name.rfind('~')
                first_end = study_name.rfind('_')
                self._name_key = study_name[0:first_end]
                self._last_name = study_name[0:last_end]
                self._first_name = study_name[last_end + 2:first_end]
                # Assign beginning and end times
                first_studies = XLTEK_Study(study_name, study_dir, False)
                self.beginning =  first_studies.creation_time
                self._end_time = first_studies._end_time
                print('%s Found.'%(self.subj))
                print('Start Time: '+str(self.beginning))
                return
        raise NameError("Could not find study with the ID : %s"% (study_ID))

    def get_studies(self):
        '''
        get_studies: Adds all studies for the subject to the studies list.
        '''
        print('Collecting Studies for %s...'%(self.subj))
        # Get all the patient's studies and store all the studies names
        study_glob = []
        for location in self.neuroworks_dirs:
            study_glob += location.glob((self._name_key + '*'))
        self._all_study_dirs |= {(os.path.split(x)[-1], x) for x in study_glob if x.is_dir()}

        # Create studies objects
        print('Creating Studies...')
        print('%d Folders To Parse...'%(len(self._all_study_dirs)))
        new_studies = self.create_studies(self._all_study_dirs, self.beginning)
        for i in range(len(new_studies)):
            new_studies[i].set_study_number(i)
            self.studies.append(new_studies[i])

        # Set the Beginning and End to their appropriate time points
        self.beginning = self.studies[0].creation_time
        print('%d Studies found for %s with %d entries.' % (len(self.studies),self.subj,self.total_data_entries))

    def create_studies(self, study_dirs, current):
        '''
        create_studies: Creates a list of study objects that are valid with in a time-frame. 
        
        Parameters
        :param study_dirs:  list        A list of path objects that are the paths to the study directories.
        :param current:     datetime    The start time point where all proceeding studies are valid.
        
        Returns
        :new_studies:   list    A list of the new studies revelent to this subject. 
        '''
        c_date = current
        to_del = []
        studies = []
        new_studies = []
        temp_studies = []

        # Get time information by creating temporary study objects
        for (name, direct) in study_dirs:
            study = XLTEK_Study(name, direct, False)
            if study.valid:
                temp_studies.append(study)
        temp_studies.sort(key=lambda study: study.creation_time)

        # Remove Studies not within the overall study time-frame
        for study in temp_studies:
            if (isinstance(study.creation_time, datetime.datetime)) and study.creation_time >= c_date:
                if len(studies) < 1 and study.study_id != self.f_study_ID:
                    continue
                if (study.creation_time - c_date) < datetime.timedelta(hours=6):
                    studies.append(study)
                    c_date = study._end_time
                else:
                    break

        # Remove Studies that overlap
        if self.studies:
            studies.insert(0,self.studies[-1])
            to_del = [0]
        if len(studies) > 1:
            for i in range(1,len(studies)):
                if studies[i-1].creation_time < studies[i].creation_time < studies[i-1]._end_time:
                    to_del.append(i)

        valid_studies = [i for j,i in enumerate(studies) if j not in to_del]

        # Build full study objects with multiprocessing
        print('Building ' + str(len(valid_studies)) + ' Studies.')
        if len(valid_studies) > 1:
            processors = len(valid_studies)-1 if len(valid_studies)<6 else 5
            with multiprocessing.Pool(processes=processors) as pool:
                study = valid_studies.pop()
                results = [pool.apply_async(XLTEK_Study, (st.study_name,st.study_dir)) for st in (valid_studies)]
                last_study = XLTEK_Study(study.study_name,study.study_dir)
                for res in results:
                    new_studies.append(res.get(timeout=120))
                new_studies.append(last_study)
            print('Sorting Studies')
            new_studies.sort(key= lambda study: study.creation_time)
        elif len(valid_studies) == 1:
            study = valid_studies.pop()
            new_studies.append(XLTEK_Study(study.study_name,study.study_dir))

        return new_studies

    def get_total_entries(self):
        '''
        get_total_entries: Gets the total number of data entries across all studies. 
        
        Returns
        :total_entries:     int     The total number of data entries across all studies. 
        '''
        total_entries = 0
        if len(self._studies) != len(self.studies):
            self._studies = list(self.studies)
        for study in self._studies:
            total_entries += study.total_entries
        return total_entries

    def get_last_ent(self):
        '''
        get_last_ent: Returns the current last data entry index for the subject.
        
        Returns
        :full_entry:    tuple   The last data entry index for the subject. (study, etc, etc-sub-entry, total) 
        '''
        total_ent = (self.total_data_entries,)
        entry = self._studies[-1].last_entry
        full_entry = (len(self._studies)-1,) + entry[:2] + total_ent
        return full_entry

    def find_entry(self, dt):
        '''
        find_entry: Returns the data entry index closest to the given datetime rounded down.
        
        Parameter 
        :param dt:      datetime    The time used to find the data entry index.
        
        Returns 
        :full_entry:    tuple       The data entry index closest to the datetime rounded down. (study, etc, etc-sub-entry, total)
        '''
        # Localize the datetime
        if dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None:
            dt = ucsf_time.localize(dt)
        
        # Determine if the requested time is within the subject's time-frame
        if self.beginning <= dt <= self.end_time:
            total_ent = 0
            # Check all studies if the time lies with in them
            for i in range(0,len(self._studies)):
                entry = self._studies[i].find_entry(dt)
                if entry:
                    full_entry = (i,) + entry[:2] +(entry[2]+total_ent,)
                    return full_entry
                else:
                    total_ent += self._studies[i]._total_entries
        return None

    def real_time_info(self):
        '''
        real_time_info: Compares the timing information of the last entry to now.

        Returns
        :diff:      timedelta    The difference between the SNC time and now.
        '''
        etc, ent, t_ent = self._studies[-1].last_entry
        entries = self._studies[-1].erds[etc].read_entries((ent,))
        now = ucsf_time.localize(datetime.datetime.now())
        diff = now - entries[0]['snc_end']
        print('SNC Time: %s'%(str(entries[0]['snc_end'])))
        print('End Time: %s'%(str(entries[0]['end_time'])))
        print('SNC/End Diff: %s'%(str(entries[0]['end_time'] - entries[0]['snc_end'])))
        print('Now: %s'%(str(now)))
        print('Now/SNC Diff: %s'%(str(diff)))
        return diff

    def get_data_range(self, start_ent, end_ent, t2d=True, buffer=None, interrupt=multiprocessing.Event()):
        '''
        get_data_range: Gets all the data between two data entry indices and either returns the data or puts it in a buffer.
        
        Parameters
        :param start_ent:   tuple       The data entry index to start getting data from.
        :param end_ent:     tuple       The last data entry index to get data from.
        :param t2d:         boolean     Return the data as either an object [False] or dictionary [True].
        :param buffer:      list        The multiprocess.Manager list append the data to.
        :param interrupt:   Event       An event that can cancel this function.

        Returns
        :data:              list        All the data in a long list if no buffer was added otherwise is None.
        '''
        data = []
        if isinstance(start_ent, tuple):
            # print('start %d end %d'%(start_ent[0],end_ent[0]))
            if start_ent[0] == end_ent[0]:
                data.extend(self._studies[start_ent[0]].get_data_range(start_ent[1:], end_ent[1:], t2d, buffer, interrupt))
            else:
                data.extend(self._studies[start_ent[0]].get_to_last(start_ent[1:], t2d, buffer, interrupt))
                for i in range(start_ent[0]+1, end_ent[0]):
                    data.extend(self._studies[i].get_to_last((0,0,0), t2d, buffer, interrupt))
                    if interrupt.is_set():
                        print('got out top')
                        break
                data.extend(self._studies[end_ent[0]].get_data_range((0,0,0), end_ent[1:], t2d, buffer, interrupt))
        if interrupt.is_set():
            print('got out top')
        return data

    def get_last_data(self, entries, t2d=True, buffer=None, interrupt=multiprocessing.Event()):
        '''
        get_last_data: Gets a specified number of data entries from the most recent data entry and either returns the
            data or puts it in a buffer.

        Parameters
        :param entries:     int         The amount of entries to get from the last entry.
        :param t2d:         boolean     The format to return the data as either a object or dictionary.
        :param buffer:      list        The multiprocess.Manager list append the data to.
        :param interrupt:   Event       An event that can cancel this function.

        Returns
        :data:              list        All the data in a long list if no buffer was added otherwise is None.
        '''
        data = []
        self.current_study = len(self._studies)
        index = 0
        study_ents = []
        # Make a list of how many entries to take from each study.
        for i in range(self.current_study-1, -1, -1):
            study_entries = self._studies[i]._total_entries
            if entries > study_entries:
                study_ents.append(study_entries)
                entries -= study_entries
            else:
                study_ents.append(entries)
                index = i
                break

        # Get data from each study
        for i in range(index, self.current_study):
            data += self._studies[i].get_last_data(study_ents.pop(-1), t2d, buffer, interrupt)
            if interrupt.is_set():
                break
        if interrupt.is_set():
            print('got out top')

        return data

    def load_buffer(self):
        '''
        load_buffer: Fill the buffer with the buffer_limit's worth of data.
        '''
        print('Filling Buffer...')
        self._buffer = self.get_last_data(self.buffer_limit)
        print('Done with filling.')

    def update_study_dirs(self):
        '''
        update_study_dirs: Checks the Neuroworks folders for new studies and adds them if valid.
        '''
        # Make a set of all patient dirs currently
        new_studies =[]
        study_glob = []
        study_dirs = set()
        for location in self.neuroworks_dirs:
            study_glob += location.glob((self._name_key + '*'))
        study_dirs |= {(os.path.split(x)[-1],x) for x in study_glob if x.is_dir()}

        # Compare the set of current dirs to ones on record
        if study_dirs != self._all_study_dirs:
            new_dirs = study_dirs - self._all_study_dirs
            if len(new_dirs) > 0:
                print('New Study: %s'%(str(datetime.datetime.now())))
                self._all_study_dirs |= new_dirs
                time.sleep(20)
                new_studies = self.create_studies(new_dirs, self.end_time)
                if len(new_studies)>0:
                    print('Added %d'%(len(new_studies)))
                    for i in range(len(new_studies)):
                        new_studies[i].set_study_number(i + len(self.studies))
                    self.studies.extend(new_studies)
                else:
                    print('Not Added')
        return new_studies

    def update_studies(self, indices):
        '''
        update_studies: Updates the etcs and erds within studies.

        Parameters
        :param indices:     list     A list of the indices of the study to update
        '''
        if type(indices) is int: indices = [indices]
        for i in indices:
            self.studies[i].get_etcs()
            self.studies[i].get_erds()

    def start_updater(self):
        '''
        start_updater: Starts the study updating in a new process.
        '''
        if not self.studies_process.is_alive():
            self.stop_events['update'].clear()
            self.studies_process = multiprocessing.Process(target=self._studies_updater, name='Update Studies')
            self.studies_process.daemon = True
            self.studies_process.start()

    def stop_updater(self):
        '''
        stop_updater: Stops the study updating process.
        '''
        if self.studies_process.is_alive():
            self.stop_events['update'].set()
            self.studies_process.join()
            self.stop_events['update'].clear()

    def start_current(self):
        '''
        start_current: Starts the current data collecting process.

        Returns
        :_buffer:       list    A multiprocessing.Managner list that the data is added to.
        '''
        if not self.current_process.is_alive():
            self.stop_events['current'].clear()
            self.current_process = multiprocessing.Process(target=self._current_process, name='Current Process',
                                                           args=(self._buffer, self.stop_events['current'], True))
            self.current_process.start()
        return self._buffer

    def stop_current(self):
        '''
        stop_current: Stops the current data collection process.
        '''
        if self.current_process.is_alive():
            self.stop_events['current'].set()
            self.current_process.join()
            self.stop_events['current'].clear()

    def start_forward(self, start_entry, spawn=False, t2d=True, buffer=None, stop_event=multiprocessing.Event()):
        '''
        start_forward: Spawns a new forward data collection process. Which takes a start data entry index and returns
            all the data from that point to the current data and continues to return new data.

        Parameters
        :param start_entry:     tuple       The data entry index to start collecting data from.
        :param spawn:           boolean     Indicates whether to make a new buffer or use the default one.
        :param t2d:             boolean     Return the data either as an object [False] or a dictionary [True]
        :param buffer:          list        A multiprocessing.Managner list to add the data to.
        :param stop_event:      Event       An event used to stop the process.

        Returns
        :buffer:        list        A multiprocessing.Managner list that the data is added to.
        :process:       object      The multiprocessing.Process object representing the process.
        :stop_event:    Event       The Event used to stop the process.
        '''
        if spawn:
            buffer = self.manager.list()
        elif buffer is None:
            buffer = self._for_buffer
        self.stop_events['forward'].append(stop_event)
        process = multiprocessing.Process(target=self._forward_process, name='Forward Process', args=(start_entry, buffer, stop_event, t2d))
        self.forward_process.append(process)
        process.start()
        return buffer, process, stop_event

    def stop_all_forward(self):
        '''
        stop_all_forward: Stops all forward data collection processes.
        '''
        self.stop_events['fill'].set()

    def shutdown(self):
        '''
        shutdown: Stops all processes.
        '''
        self.stop_all_forward()
        self.stop_current()
        self.stop_updater()
        self.manager.shutdown()

    def _studies_updater(self):
        '''
        _studies_updater: The loop function that the studies updater process will run.
        :return:
        '''
        print('Running Updater')
        # Continuously check if there is are new studies and add them if there are.
        while not self.stop_events['update'].is_set():
            self.update_study_dirs()

    def _buffer_process(self, buffer, local_entries=None, t2d=True, stop_event=multiprocessing.Event(), finfout=True):
        '''
        _buffer_process: The loop function that adds new data to a buffer.

        Parameters:
        :param buffer:          list        The multiprocessing.Manager list to add the data to.
        :param local_entries:   int         The total entry number to start collecting from.
        :param t2d:             boolean     Return the data either as an object [False] or a dictionary [True]
        :param stop_event:      Event       The Event the used to stop the loop.
        :param finfout:         boolean     Controls if the buffer is first in first out [True] or waits until there is space [False].
        '''
        # print('buffering')
        # Set the entry to start collecting from.
        if local_entries is None:
            local_entries = self.total_data_entries
        # Continuously check if there is a difference in the number of entries and if there is collect the new data
        while not (self.stop_events['fill'].is_set() or stop_event.is_set()):
            new_total = self.total_data_entries
            new_entries = new_total - local_entries
            if new_entries > 0:
                if finfout:
                    if new_entries < self.buffer_limit:
                        if len(buffer) + new_entries > self.buffer_limit:
                            del buffer[0:new_entries]
                        self.get_last_data(new_entries, t2d, buffer)
                    else:
                        buffer.clear()
                        self.get_last_data(self.buffer_limit, t2d, buffer)
                    local_entries = new_total
                else:
                    if new_entries < self.buffer_limit:
                        if len(buffer) + new_entries > self.buffer_limit:
                            continue
                    self.get_last_data(new_entries, t2d, buffer)
                    local_entries = new_total

    def _forward_process(self, start_ent, buffer, stop_event, t2d=True):
        '''
        _forward_process: Gets all data from a data entry index to the current time and adds any new data collected.
            Also, the buffer will fill up to the limit then pause until the data is removed from the buffer.

        Parameters
        :param start_ent:       tuple       The data entry index of the first data entry to grab.
        :param buffer:          list        The multiprocessing.Manager list that the data is added to.
        :param stop_event:      Event       The Event used to stop the loop.
        :param t2d:             boolean     Return the data either as an object [False] or a dictionary [True]
        '''
        end_ent = self.get_last_ent()
        self.get_data_range(start_ent, end_ent, t2d, buffer, stop_event)
        print('A forward process has returned all previous data and will now return incoming data.')
        self._buffer_process(buffer, end_ent[3], t2d, stop_event, False)
        
    def _current_process(self, buffer, stop_event, t2d=True):
        '''
        _current_process: The loop function adds the most recent data up to the buffer limit and updates the data. When
            the buffer limit is reached the new data will override the oldest data.

        Parameters
        :param buffer:          list        The multiprocessing.Manager list that the data is added to.
        :param stop_event:      Event       The Event used to stop the loop.
        :param t2d:             boolean     Return the data either as an object [False] or a dictionary [True]
        '''
        local_entries = self.total_data_entries
        self.get_last_data(self.buffer_limit, t2d, buffer)
        print('Buffer has been Filled')
        self._buffer_process(buffer, local_entries, t2d, stop_event)

    def proxy(self, name):
        '''
        proxy: Changes an attribute of this object to a proxy. Meaning it references a different dictionary instead of
            its __dict__.

        Parameter
        :param name:    string      The name of the attribute to change.
        '''
        def get(self):
            return self.proxy_dict[name]
        def set(self, value):
            self.proxy_dict[name] = value
        setattr(XLTEK_Reader, name, property(get, set))

    def proxies(self, exclude=('proxy_dict',)):
        '''
        proxies: Changes all attributes of this object to proxies except ones specified. Meaning the all chosen attributes
            will references a different dictionary instead of this object's __dict__.

        Parameters
        :param exclude:     tuple   A tuple of the names of all the attributes that will not have proxies.
        '''
        remove = []
        for key in self.__dict__:
            if key not in exclude:
                self.proxy_dict[key] = self.__dict__[key]
                remove.append(key)
        for key in remove:
            self.__dict__.pop(key)
            self.proxy(key)

    @property
    def end_time(self):
        '''
        end_time: An property function that return the end timestamp of last data of the last study.

        Returns
        :return:    datetime     The end timestamp of last data of the last study.
        '''
        if len(self._studies) != len(self.studies):
            self._studies = list(self.studies)
        self._end_time = self._studies[-1].end_time
        return self._end_time

    @property
    def total_data_entries(self):
        '''
        total_data_entries: An property function that return the total number of entries across all studies.

        Returns
        :return:    int      The total number of entries across all studies.
        '''
        return self.get_total_entries()


class XLTEK_Study:
    """
    XLTEK_Study: An object that contains a study created by the XLTEK system.

    Required Modules: os, pathlib, struct, re, multiprocessing, time, datetime, pytz, xlrdlib2, numpy
    Required Classes: StcFile, SncFile, EtcFile, ErdFile
    Methods: validate, get_time_info, set_study_number, build, get_etcs, get_erds, create_etcs, create_erds, find_entry,
        get_data_range, get_to_last, get_last_data, update_entries, snc_sample_stamp

    Class Attributes
    none

    Object Parameters & Attributes
    Parameters
    :param study_name:  string      The full XLTEK name of the study.
    :param study_dir:   string      The directory of the study.
    :full_init:         boolean     Determines if the etcs and erds are loaded.

    Attributes
    :header:            object      The XLTEK header of the stc file.
    :study_name:        string      The name of the XLTEK study.
    :study_number:      int         The number of the study for the subject.
    :study_dir:         Path        The path of the the study's directory.
    :complete:          boolean     A flag for determining if the study is complete or still running.
    :valid:             boolean     A flag for determining if the study has all the necessary components.

    :creation_time:     datetime    The datetime of the start of the study.
    :end_time:          datetime    The datetime of the end of the study.
    :sample_offset:     int         The number of empty samples from the beginning before data was recorded.
    :start_sample:      int         The sample number of the first sample.
    :end_sample:        int         The sample number of the last sample.
    :sample_rate:       int         The sampling rate of the study.
    :last_entry:        tuple       The last data entry index of the study.
    :total_entries:     int         The total number of entries for the study with updating.
    :_total_entries:    int         The total number of entries for the study without updating.

    :etc_paths:         set         The paths of all the etcs of the study.
    :erd_paths:         set         The paths of all the erds of the study.

    :stc:               Stcfile     The object that maps the etc and erd files.
    :snc:               Sncfile     The object that maps samples to wall clock time.
    :epo:               Path        The path of the epo file. (unused right now)
    :etcs:              list        A list of EtcFile that contains all of the sample information about the data
    :erds:              list        A list of ErdFile that holds all the data.
    """
    def __repr__(self):
        return repr((self.creation_time, self.end_time, self.start_sample, self.end_sample))

    def __init__(self, study_name, study_dir, full_init=True):
        # Identifiers #
        self.header = None
        self.study_name = study_name
        self.study_id = study_name.split('_')[-1]
        self.study_number = 0
        if isinstance(study_dir, pathlib.Path):
            self.study_dir = study_dir
        else:
            self.study_dir = pathlib.Path(study_dir)
        self.complete = False

        # Time Information #
        self.creation_time = None
        self._end_time = None
        self.sample_offset = 0
        self.start_sample = None
        self.end_sample = None
        self.sample_rate = 0

        # File Paths #
        self.etc_paths = set()
        self.erd_paths = set()

        # Study Objects #
        self.stc = None
        self.snc = None
        self.epo = None
        self.etcs = []
        self.erds = []

        ## Initialization Methods ##
        self.valid = self.validate()

        if self.valid:
            self.get_time_info()
            if full_init:
                self.build()


    def validate(self):
        '''
        validate: Checks if the STC file and the SNC file are present and will wait 10 seconds before returning False.

        Returns
        :flag:      boolean     Tell whether the files are present or not.
        '''
        stc_loc = pathlib.Path(self.study_dir, self.study_name + '.stc')
        snc_loc = pathlib.Path(self.study_dir, self.study_name + '.snc')
        # Wait for the file to be generated if it is not there.
        flag1 = file_wait(stc_loc, 10)
        flag2 = file_wait(snc_loc, 10)
        flag = flag1 and flag2
        if flag:
            self.stc = StcFile(stc_loc)
            self.snc = SncFile(snc_loc)
            # If there is no data then the flag is set to False.
            if len(self.stc.stc_data) <= 0:
                flag = False
        return flag

    def get_time_info(self):
        '''
        get_time_info: Reads the STC and SNC files to find the beginning and end of study.
        '''
        self.header = self.stc._generic_hdr
        self.creation_time = self.header.creation_time
        f_erd = ErdFile(pathlib.Path(self.study_dir, self.stc.stc_data[0].segment_name+'.erd'))
        self.sample_rate = f_erd.header.sample_freq
        self.sample_offset = self.stc.stc_data[0].start_stamp
        self.start_sample = self.stc.stc_data[0].start_stamp
        self.end_sample = self.stc.stc_data[-1].end_stamp
        self._end_time = self.header.creation_time + datetime.timedelta(seconds=(self.end_sample)/self.sample_rate)
        if self.creation_time > self.snc.snc_data[0].sample_time:
            self.creation_time = self.snc.snc_data[0].sample_time
        if ucsf_time.localize(datetime.datetime.now()) - self._end_time > datetime.timedelta(hours=1):
            self._end_time = self.snc.snc_data[-1].sample_time
            self.complete = True

    def build(self):
        '''
        build: Gets the ETC files, ERD files, sets the current ETC, and the total number of entries.
        '''
        self.get_etcs()
        self.get_erds()

        self.current_etc = len(self.etcs)
        self.etc_total = len(self.etcs)
        self.last_etc_size = len(self.etcs[self.etc_total - 1].etc_data)
        self._subtotal_entries = 0
        for i in range(0, self.etc_total - 1):
            self._subtotal_entries += len(self.etcs[i].etc_data)

    def set_study_number(self, number):
        '''
        set_study_number: Sets the study number for all etc and erd objects and assigns their index number.

        Parameters
        :param number:      int     The number of the study to assign this study to.
        '''
        self.study_number = number
        for i in range(0, len(self.etcs)):
            self.etcs[i].study_number = self.study_number
            self.etcs[i].etc_number = i
        for i in range(0, len(self.erds)):
            self.erds[i].study_number = self.study_number
            self.erds[i].erd_number = i

    def get_etcs(self):
        '''
        get_etcs: Gets ETCs and checks if there are new ones.
        '''
        # Get the file names
        etc_paths = {pathlib.Path(self.study_dir, etc.segment_name+'.etc') for etc in self.stc.stc_data}

        # Compare the set of current file names to ones on record
        #print('etc: %d stc: %d self: %d'%(len(etc_paths),len(self.stc.stc_data),len(self.etc_paths)))
        if etc_paths != self.etc_paths:
            new_etcs = etc_paths - self.etc_paths
            self.etc_paths |= new_etcs
            self.etcs.extend(self.create_etcs(new_etcs))
            self.etcs.sort(key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', str(s._fpath))])
            for i in range(0,len(self.etcs)):
                self.etcs[i].etc_number = i

    def get_erds(self):
        '''
        get_erds: Gets the ERDs and checks if there are new ones.
        '''
        # Get the file names
        erd_paths = {pathlib.Path(self.study_dir, etc.segment_name+'.erd') for etc in self.stc.stc_data}
        # Compare the set of current file names to ones on record
        if erd_paths != self.erd_paths:
            new_erds = erd_paths - self.erd_paths
            self.erd_paths |= new_erds
            self.erds.extend(self.create_erds(new_erds))
            self.erds.sort(key=lambda s: [int(t) if t.isdigit() else t.lower() for t in re.split('(\d+)', str(s._fpath))])
            for i in range(0,len(self.erds)):
                self.erds[i].erd_number = i

    def create_etcs(self, etc_fpaths):
        '''
        create_etcs: Makes EtcFile objects.

        Parameters
        :param etc_fpaths:      list    A list of paths to ETC files.

        Returns
        :etcs:                  list    A list of the new EtcFiles.
        '''
        etcs = []
        for etc_fpath in etc_fpaths:
            new_etc = EtcFile(etc_fpath, self.snc, self.study_number)
            etcs.append(new_etc)
        return etcs

    def create_erds(self, erd_fpaths):
        '''
        create_erds: Makes ErdFile objects.

        Parameters
        :param erd_fpaths:      list    A list of paths to ERD files.

        Returns
        :erds:                  list    A list of the new ErdFiles.
        '''
        erds = []
        for erd_fpath in erd_fpaths:
            etc = next((x for x in self.etcs if str(x._fpath)[:-4] == str(erd_fpath)[:-4]), None)
            new_erd = ErdFile(erd_fpath, self.snc, self.creation_time, self.sample_offset, etc, self.study_number)
            erds.append(new_erd)
        return erds

    def update_entries(self):
        '''
        update_entries: Updates the etc files, erd files, and the number of total entries.
        '''
        # Update the etcs and erds
        if not self.complete:
            self.get_etcs()
            self.get_erds()
        self.current_etc = len(self.etcs)

        # Update the subtotal number of entries
        if self.current_etc != self.etc_total:
            self._subtotal_entries += len(self.etcs[self.etc_total - 1].etc_data)
            self.etc_total = self.current_etc
        self.last_etc_size = len(self.etcs[self.etc_total - 1].etc_data)

    def find_entry(self, dt):
        '''
        find_entry: Finds the nearest data entry index to a given time.

        Parameters
        :param dt:  datetime    The datetime to find the data entry index for.

        Returns
        :entry:     tuple       The data entry index closest to the given time.
        '''
        # Set the timezone if there is none
        if dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None:
            dt = ucsf_time.localize(dt)
        # Check if the time is within the study
        if self.creation_time <= dt <= self.end_time:
            distance = 0
            self.update_entries()
            # For each etc file check if the time is within it.
            for i in range(0, self.etc_total):
                etc_start = self.stc.stc_data[i].start_stamp - self.sample_offset
                etc_end = etc_start + self.stc.stc_data[i].sample_span
                start = datetime.timedelta(seconds=(etc_start / self.sample_rate))
                end = datetime.timedelta(seconds=(etc_end / self.sample_rate))
                data_start = self.creation_time + start
                if data_start <= dt <= self.creation_time + end:
                    # For each entry check if the time is within it.
                    for j in range(0, len(self.etcs[i].etc_data)):
                        entry_begin = self.etcs[i].etc_data[j].sample_num
                        entry_last = entry_begin + self.etcs[i].etc_data[j].sample_span
                        begin = datetime.timedelta(seconds=(entry_begin / self.sample_rate))
                        last = datetime.timedelta(seconds=(entry_last) / self.sample_rate)
                        if data_start + begin <= dt < data_start + last:
                            return i, j, (distance + j)
                else:
                    distance += len(self.etcs[i].etc_data)
        return None

    def find_sync_sample(self, time_stamp):
        '''
        find_sync_sample: Finds the nearest sample to a given time.

        Parameters
        :param time_stamp:  datetime    The datetime to find the sample for.

        Returns
        :entry:             int         The sample closest to the given time. Returns -1 if it is before the study.
        '''
        if self.snc.snc_data[0].sample_time > time_stamp:
            return -1
        elif self.snc.snc_data[-1].sample_time < time_stamp:
            return self.snc.snc_data[-1].sample_stamp
        else:
            for i in range(0, len(self.snc.snc_data) - 1):
                if self.snc.snc_data[i].sample_time < time_stamp < self.snc.snc_data[i + 1].sample_time:
                    return self.snc.snc_data[i].sample_stamp

    def get_data_range(self, start_ent, end_ent, t2d=True, buffer=None, interrupt=multiprocessing.Event()):
        '''
        get_data_range: Gets the data between two data entry indices and either returns it or sends it to a buffer.

        Parameter
        :param start_ent:       tuple       The data entry index to start getting data from. (etc, entry, total entry)
        :param end_ent:         tuple       The last data entry index to get data from. (etc, entry, total entry)
        :param t2d:             boolean     Return the data as either an object [False] or dictionary [True].
        :param buffer:          list        The multiprocess.Manager list append the data to.
        :param interrupt:       Event       An event that can cancel this function.

        Returns
        :data:                  list        All the data in a long list if no buffer was added otherwise is None.
        '''
        data = []
        if isinstance(start_ent, tuple):
            s_etc, s_entry, s_t_ent = start_ent
            e_etc, e_entry, e_t_ent = end_ent
            if s_etc == e_etc:
                data = self.erds[s_etc].read_entries(range(s_entry, e_entry), t2d, buffer, interrupt)
            else:
                data += self.erds[s_etc].read_entries(range(s_entry, len(self.etcs[s_etc].etc_data)), t2d, buffer, interrupt)
                for i in range(s_etc + 1, e_etc):
                    data += self.erds[i].read_entries(range(0, len(self.etcs[i].etc_data)), t2d, buffer, interrupt)
                    if interrupt.is_set():
                        print('escaping middle')
                        break
                data += self.erds[e_etc].read_entries(range(0, e_entry), t2d, buffer, interrupt)
        return data

    def get_to_last(self, start_ent, t2d=True, buffer=None, interrupt=multiprocessing.Event()):
        '''
        get_to_last: Gets the from a data entry index to the end and either returns it or sends it to a buffer.

        Parameter
        :param start_ent:       tuple       The data entry index to start getting data from. (etc, entry, total entry)
        :param t2d:             boolean     Return the data as either an object [False] or dictionary [True].
        :param buffer:          list        The multiprocess.Manager list append the data to.
        :param interrupt:       Event       An event that can cancel this function.

        Returns
        :data:                  list        All the data in a long list if no buffer was added otherwise is None.
        '''
        data = []
        if isinstance(start_ent, tuple):
            etc, entry, t_ent = start_ent
            data = self.erds[etc].read_entries(range(entry, len(self.etcs[etc].etc_data)), t2d, buffer)
            for i in range(etc+1,self.etc_total):
                if len(self.etcs[i].etc_data) > 0:
                    data.extend(self.erds[i].read_entries(range(0,len(self.etcs[i].etc_data)), t2d, buffer))
                if interrupt.is_set():
                    print('escaping middle')
                    break
        return data

    def get_last_data(self, entries=1, t2d=True, buffer=None, interrupt=multiprocessing.Event()):
        '''
        get_last_data: Gets a specified number of data entries from the most recent data entry and either returns the
            data or puts it in a buffer.

        Parameter
        :param entries:     int         The amount of entries to get from the last entry.
        :param t2d:         boolean     Return the data as either an object [False] or dictionary [True].
        :param buffer:      list        The multiprocess.Manager list append the data to.
        :param interrupt:   Event       An event that can cancel this function.

        Returns
        :data:              list        All the data in a long list if no buffer was added otherwise is None.
        '''
        data = []
        if entries <= self.last_etc_size:
            data += self.erds[self.etc_total-1].read_entries(range(self.last_etc_size-entries, self.last_etc_size), t2d, buffer)
        else:
            index = 0
            ranges = [range(0,self.last_etc_size)]
            entries -= self.last_etc_size
            for i in range(self.etc_total-2, -1, -1):
                if entries > len(self.etcs[i].etc_data):
                    ranges.append(range(0, len(self.etcs[i].etc_data)))
                    entries -= len(self.etcs[i].etc_data)
                else:
                    ranges.append(range(len(self.etcs[i].etc_data)-entries, len(self.etcs[i].etc_data)))
                    index = i
                    break

            if buffer is not None:
                for i in range(index, self.etc_total):
                    data += self.erds[i].read_entries(ranges.pop(-1), t2d, buffer)
                    if interrupt.is_set():
                        print('escaping middle')
                        break
            else:
                with NormalPool(processes=4) as pool:
                    results = []
                    first = ranges.pop(-1)
                    for i in range(index, self.etc_total):
                        results.append(pool.apply_async(read_erd, (self.erds[i],ranges.pop(-1),t2d)))
                        if len(ranges) < 2:
                            break
                    data += self.erds[self.etc_total - 1].read_entries(first, t2d)
                    for res in results:
                        data += res.get(timeout=600)
                        if interrupt.is_set():
                            print('escaping middle')
                            break
        if interrupt.is_set():
            print('got out middle')
        return data

    @property
    def last_entry(self):
        '''
        last_entry: An property function that return the last data entry index of the study.

        Returns
        :return:    tuple      The last entry index of the study.
        '''
        self.update_entries()
        return (self.etc_total-1, self.last_etc_size-1, self._total_entries)

    @property
    def end_time(self):
        '''
        end_time: An property function that return the end timestamp of last data in the study.

        Returns
        :return:    datetime     The end timestamp of last data in the study.
        '''
        try:
            etc, ent, t_ent = self.last_entry
            entries = self.erds[etc].read_entries((ent,))
            self._end_time = entries[0]['snc_end']
        finally:
            return self._end_time

    @property
    def total_entries(self):
        '''
        total_entries: An property function that return the total number of entries for the study with updating.

        Returns
        :return:    int      The total number of entries for the study.
        '''
        self.update_entries()
        return (self._subtotal_entries + self.last_etc_size)

    @property
    def _total_entries(self):
        '''
        total_data_entries: An property function that return the total number of entries for the study without updating.

        Returns
        :return:    int      The total number of entries for the study.
        '''
        return (self._subtotal_entries + self.last_etc_size)


class XLRecListFileReader(object):
    def __init__(self, filepath, data_offset, schema, create_rec):
        self._fpath = filepath
        self._start_offset = data_offset
        self.file_schema = schema
        self._pos = self._start_offset

        self.new_rec = create_rec

        return

    def __iter__(self):
        return self

    def __next__(self):
        pos = self._pos

        try:
            with open(self._fpath, 'rb') as fobj:
                fobj.seek(pos)
                rec = self.new_rec(fobj, self.file_schema)
                pos = fobj.tell()
        except:
            raise StopIteration()

        self._pos = pos
        return rec

    def next(self):
        return self.__next__()

    def reset(self, data_offset=None):
        if data_offset is None:
            data_offset = self._start_offset
        self._pos = data_offset

    def list_gen(self):
        pos = self._pos
        reading = True
        recs = []
        with open(self._fpath, 'rb') as fobj:
            while reading:
                try:
                    fobj.seek(pos)
                    recs.append(self.new_rec(fobj, self.file_schema))
                    pos = fobj.tell()
                except:
                    reading = False
        self._pos = pos
        return recs

class XLTEK_GenericHeader(object):
    HDRSIZE = 352

    def __init__(self, fobj):
        self._fobj = fobj
        pos = self._fobj.tell()

        try:
            self.parse()
            self.complete = True
        except struct.error as e:
            self.complete = False
        self._fobj.seek(pos)
        return

    def parse(self):
        self._fobj.seek(0)  # should be at start of file
        self.guid = self._fobj.read(16)
        schema_bytes = self._fobj.read(4)
        if struct.unpack('I', schema_bytes) == 0:
            self.file_schema = struct.unpack('I', schema_bytes)
        else:
            (self.file_schema, self.base_schema) = struct.unpack('HH', schema_bytes)

        creation_time = datetime.datetime.fromtimestamp(float(struct.unpack('I', self._fobj.read(4))[0]))
        self.creation_time = pytz.timezone('America/Los_Angeles').localize(creation_time)
        if self.file_schema == 0:
            self.product_version = struct.unpack('II', self._fobj.read(8))
        else:
            self._patient_id, = struct.unpack('I', self._fobj.read(4))
            self._study_id, = struct.unpack('I', self._fobj.read(4))

        self.patient_name = (self._fobj.read(80),
                             self._fobj.read(80),
                             self._fobj.read(80))
        self.patient_id = self._fobj.read(80)

        self.HDRSIZE = self._fobj.tell()
        return

    def __getstate__(self):
        state = self.__dict__.copy()
        fobj = state['_fobj']
        state['_fobj'] = (fobj.name, fobj.mode)
        return state

    def __setstate__(self, state):
        name, mode = state['_fobj']
        state['_fobj'] = open(name, mode)
        state['_fobj'].close()
        self.__dict__.update(state)

class StcFile(object):
    def __init__(self,filename):
        """
        :param filename:
        :return:
        """
        if isinstance(filename, pathlib.Path):
            self._fpath = filename
        else:
            self._fpath = pathlib.Path(filename)

        file_wait(self._fpath)
        with open(self._fpath, 'rb') as fobj:
            self._generic_hdr = XLTEK_GenericHeader(fobj)
            self.complete = self._generic_hdr.complete
            if self.complete:
                self.stc_header = StcHeader(fobj)

        if self.complete:
            self._data_offset = self.stc_header.STC_DATA_OFFSET
            self._reader = XLRecListFileReader(self._fpath, self._data_offset, self._generic_hdr.file_schema, new_stc_entry)

        self._stc_data = None
        return

    def __iter__(self):
        return iter(self.stc_data)

    def build(self):
        file_wait(self._fpath)
        with open(self._fpath, 'rb') as fobj:
            self._generic_hdr = XLTEK_GenericHeader(fobj)
            self.complete = self._generic_hdr.complete
            if self.complete:
                self.stc_header = StcHeader(fobj)

        if self.complete:
            self._data_offset = self.stc_header.STC_DATA_OFFSET
            self._reader = XLRecListFileReader(self._fpath, self._data_offset, self._generic_hdr.file_schema, new_stc_entry)

    def read_data(self):
        if not self._stc_data:
            self._stc_data = self._reader.list_gen()
        else:
            self._stc_data.extend(self._reader.list_gen())
        return self._stc_data

    @property
    def stc_data(self):
        return self.read_data()

class StcHeader(object):
    STC_DATA_OFFSET = 408

    def __init__(self, fobj):
        self._fobj = fobj

        self._fobj.seek(XLTEK_GenericHeader.HDRSIZE)
        (self.next_segment,) = struct.unpack('I',self._fobj.read(4))  # these don't appear accurate
        (self.final,) = struct.unpack('I',self._fobj.read(4))         # ibid.
        self.padding = list(struct.unpack('12I',self._fobj.read(48))) # ibid.

        self._data_offset = self._fobj.tell()

        #assert(self._data_offset == self.STC_DATA_OFFSET)
        return

    def __getstate__(self):
        state = self.__dict__.copy()
        fobj = state['_fobj']
        state['_fobj'] = (fobj.name, fobj.mode)
        return state

    def __setstate__(self, state):
        name, mode = state['_fobj']
        state['_fobj'] = open(name, mode)
        state['_fobj'].close()
        self.__dict__.update(state)

class StcEntry(object):
    def __init__(self, fobj, fschema):
        self.segment_name = fobj.read(256).decode('ascii').strip('\x00')
        (self.start_stamp,) = struct.unpack('I', fobj.read(4))
        (self.end_stamp,) = struct.unpack('I', fobj.read(4))
        (self.sample_num,) = struct.unpack('I', fobj.read(4))
        (self.sample_span,) = struct.unpack('I', fobj.read(4))

        return

class SncFile(object):
    def __init__(self, filename):
        if isinstance(filename, pathlib.Path):
            self._fpath = filename
        else:
            self._fpath = pathlib.Path(filename)

        # headers
        with open(self._fpath, 'rb') as fobj:
            self._generic_hdr = XLTEK_GenericHeader(fobj)
            self._data_offset = self._generic_hdr.HDRSIZE

        self._reader = XLRecListFileReader(self._fpath, self._data_offset, self._generic_hdr.file_schema, new_snc_entry)

        self._snc_data = None
        return

    def __iter__(self):
        return iter(self.snc_data)

    def read_data(self):
        if not self._snc_data:
            self._snc_data = self._reader.list_gen()
        else:
            self._snc_data.extend(self._reader.list_gen())
        return self._snc_data

    @property
    def snc_data(self):
        return self.read_data()

class SncEntry(object):
    def __init__(self, fobj, fschema):
        (self.sample_stamp,) = struct.unpack('I', fobj.read(4))
        self.sample_time = filetime_to_datetime(struct.unpack('II', fobj.read(8)))

        return

class EtcFile(object):
    def __init__(self, filename, snc_file, study=0):
        if isinstance(filename, pathlib.Path):
            self._fpath = filename
        else:
            self._fpath = pathlib.Path(filename)
        self.study_number = study
        self.etc_number = 0

        self._reader = None
        self.snc = snc_file
        file_wait(self._fpath)
        with open(self._fpath, 'rb') as fobj:
            self._generic_hdr = XLTEK_GenericHeader(fobj)
        self.valid = self._generic_hdr.complete

        if self.valid:
            self._data_offset = self._generic_hdr.HDRSIZE
            self._reader = XLRecListFileReader(self._fpath, self._data_offset, self._generic_hdr.file_schema, new_etc_entry)

        self._etc_data = None
        return

    def __repr__(self):
        return repr((self._fpath))

    def __iter__(self):
        return iter(self.etc_data)

    def build(self):
        with open(self._fpath, 'rb') as fobj:
            self._generic_hdr = XLTEK_GenericHeader(fobj)
        self.valid = self._generic_hdr.complete

        if self.valid:
            self._data_offset = self._generic_hdr.HDRSIZE
            self._reader = XLRecListFileReader(self._fpath, self._data_offset, self._generic_hdr.file_schema, new_etc_entry)

    def add_snc_times(self, etcs):
        snc_len = len(self.snc.snc_data)
        for i in range(0,len(etcs)):
            start = etcs[i].samplestamp
            if snc_len > 0:
                if start < self.snc._snc_data[0].sample_stamp:
                    etcs[i].snc_sample = 1
                    etcs[i].snc_time = self.snc._generic_hdr.creation_time
                elif start >= self.snc._snc_data[-1].sample_stamp:
                    etcs[i].snc_sample = self.snc._snc_data[-1].sample_stamp
                    etcs[i].snc_time = self.snc._snc_data[-1].sample_time
                else:
                    for j in range(1, snc_len):
                        if start < self.snc._snc_data[j].sample_stamp:
                            etcs[i].snc_sample = self.snc._snc_data[j-1].sample_stamp
                            etcs[i].snc_time = self.snc._snc_data[j-1].sample_time
                            break
            else:
                etcs[i].snc_sample = 1
                etcs[i].snc_time = self.snc._generic_hdr.creation_time
        return etcs

    def read_data(self):
        if self._reader is None:
            self._etc_data = []
        elif not self._etc_data:
            temp = self._reader.list_gen()
            self._etc_data = self.add_snc_times(temp)
        else:
            temp = self._reader.list_gen()
            self._etc_data.extend(self.add_snc_times(temp))
        return self._etc_data

    @property
    def etc_data(self):
        return self.read_data()

class EtcEntry(object):
    def __init__(self, fobj, fschema):
        (self.offset,) = struct.unpack('I', fobj.read(4))
        if fschema == 2:
            (self.timestamp,) = struct.unpack('I', fobj.read(4))
        elif fschema == 3:
            (self.samplestamp,) = struct.unpack('I', fobj.read(4))
        else:
            raise NotImplementedError
        (self.sample_num,) = struct.unpack('I', fobj.read(4))
        self.sample_span = struct.unpack('HBB', fobj.read(4))[0]
        self.snc_sample = -1
        self.snc_time = None

        return

class ErdFile(object):
    def __init__(self, filename, snc_file=None, study_start=None, study_offset=None, etc=None, study=0):
        file_wait(filename)
        self._fobj = open(filename, 'rb')
        self._generic_hdr = XLTEK_GenericHeader(self._fobj)
        self.header = ErdHeader(self._fobj, study_start, study_offset, generic_hdr=self._generic_hdr)
        self.study_number = study
        self.erd_number = 0
        self.complete = self._generic_hdr.complete and self.header.complete

        if isinstance(filename, pathlib.Path):
            self._fpath = filename
        else:
            self._fpath = pathlib.Path(filename)
        if etc:
            self.etc_path = etc._fpath
            self.etc = etc
        else:
            self.etc_path = '{:s}.etc'.format('.'.join(self._fobj.name.split('.')[:-1]))
            self.etc = EtcFile(self.etc_path, snc_file)

        if self.complete:
            self._cf_arr = self._lookup_conversion_factor(self.header.num_channels,
                                                          self.header.headbox_type,
                                                          self.header.headbox_sw_version,
                                                          self.header.discardbits)
            self._erd_data = []

            self._pos = None
            self._suspend()

        return

    def __getstate__(self):
        state = self.__dict__.copy()
        fobj = state['_fobj']
        state['_fobj'] = (fobj.name, fobj.mode)
        return state

    def __setstate__(self, state):
        name, mode = state['_fobj']
        state['_fobj'] = open(name, mode)
        state['_fobj'].close()
        self.__dict__.update(state)

    def __repr__(self):
        return repr((self._fpath))

    def _lookup_conversion_factor(self, chans, hb_type, hb_sw_vers, discardbits):
        """Returns an array of size (1, num_channels) where each value is the conversion factor to go from stored values to physical (engineering) units"""

        # Conversion factor constants:
        # In the documentation it says to bit shift the data by m_discardbits (2**discardbits) but that is done in the
        # C code function that part is skipped.
        _FAC00 = (8711.0 / (2 ** 21 - 0.5))  # * (2**discardbits)
        _FAC01 = ((5000000.0 / (2 ** 10 - 0.5)) / 2 ** 6)  # * (2**discardbits)
        _FAC02 = ((1.0 / 2.0) ** 6)  # * (2**discardbits)
        _FAC03 = ((8711.0 / (2 ** 21 - 0.5)) / (159.8 / 249.5))  # * (2**discardbits)
        _FAC04 = ((10000000.0 / (2 ** 10 - 0.5)) / 2 ** 6)  # * (2**discardbits)
        _FAC05 = ((20000000.0 / 65536.0) / 2 ** 6)  # * (2**discardbits)
        _FAC06 = ((10800000.0 / 65536.0) / 2 ** 6)  # * (2**discardbits)
        _FAC07 = ((10000000.0 / 65536.0) / 2 ** 6)  # * (2**discardbits)

        cf_arr = [None for ii in range(len(hb_type))]  # KLUDGE (unclear if more than one hb will ever be connected) ... this code is ultimately not setup for that case
        for hb_idx, (hb, sw) in enumerate(zip(hb_type, hb_sw_vers)):
            # EEG32 or EEG32U
            if (hb == 1) or (hb == 19):
                cf_arr[hb_idx] = np.full(32, _FAC00)
            # EEG128
            elif hb == 3:
                cf_arr[hb_idx] = np.full(128, _FAC00)
            # EMU36
            elif hb == 6:
                cf_arr[hb_idx] = np.concatenate((np.full(32, _FAC00), np.full(4, _FAC01)))
            # AMB28
            elif hb == 4:
                cf_arr[hb_idx] = np.concatenate((np.full(24, _FAC00), np.full(4, _FAC01)))
            # MOBEE32
            elif hb == 9:
                cf_arr[hb_idx] = np.concatenate((np.full(33, _FAC00), np.full(2, _FAC02)))
            # MOBEE24
            elif hb == 8:
                cf_arr[hb_idx] = np.concatenate((np.full(25, _FAC00), np.full(2, _FAC02)))
            # HYPPO
            elif hb == 5:
                if float(sw) < 3.4:
                    _FACXX = _FAC04
                else:
                    _FACXX = _FAC05
                cf_arr[hb_idx] = np.concatenate((np.full(26, _FAC00), np.full(6, _FAC03), np.full(8, _FACXX), np.full(2, _FAC02)))
            # Connex
            elif hb == 14:
                cf_arr[hb_idx] = np.concatenate((np.full(38, _FAC00), np.full(10, _FAC06), np.full(2, _FAC02)))
            # Trex
            elif hb == 15:
                cf_arr[hb_idx] = np.concatenate((np.full(24, _FAC00), np.full(4, _FAC00), np.full(4, _FAC07), np.full(2, _FAC02)))
            # EMU40
            elif hb == 17:
                cf_arr[hb_idx] = np.concatenate((np.full(40, _FAC00), np.full(4, _FAC06), np.full(2, _FAC02)))
            # NeuroLink IP
            elif hb == 21:
                cf_arr[hb_idx] = np.concatenate((np.full(128, _FAC00), np.full(2, _FAC02), np.full(126, _FAC00)))
            # Netlink
            elif hb == 22:
                cf_arr[hb_idx] = np.concatenate((np.full(32, _FAC00), np.full(8, _FAC06), np.full(2, _FAC02), np.full(1, _FAC06)))
            # Traveler
            elif hb == 23:
                cf_arr[hb_idx] = np.concatenate((np.full(32, _FAC00), np.full(4, _FAC06), np.full(2, _FAC02), np.full(1, _FAC06)))
            # Quantum
            elif hb == 20:

                cf_arr[hb_idx] = np.concatenate((np.full(chans-20, _FAC00), np.full(17, _FAC06),np.full(3, _FAC02)))

        return cf_arr

    def _suspend(self):
        if not self._fobj.closed:
            self._pos = self._fobj.tell()
            self._fobj.close()
        return

    def _resume(self):
        if self._fobj.closed:
            self._fobj = open(self._fobj.name)
            self._fobj.seek(self._pos)
        return

    def toc2dict(self, toc):
        entry = {}
        entry['data'] = np.multiply(toc.data,toc._cf_arr[0])
        entry['channels'] = toc.num_channels
        entry['sample_rate'] = toc.sample_rate
        entry['study_start'] = toc.study_start
        entry['study_offset'] = toc.study_offset
        entry['start_sample'] = toc.toc_samplestamp
        entry['end_sample'] = entry['start_sample'] + toc.toc_nsamples
        entry['snc_sample'] = toc.snc_sample
        entry['snc_time'] = toc.snc_time
        entry['sample_span'] = toc.toc_nsamples
        entry['start_time'] = toc.study_start + datetime.timedelta(seconds=((entry['start_sample']) / toc.sample_rate))
        entry['end_time'] = toc.study_start + datetime.timedelta(seconds=((entry['end_sample']) / toc.sample_rate))
        entry['snc_start'] = toc.snc_time + datetime.timedelta(seconds=(entry['start_sample']-toc.snc_sample)/toc.sample_rate)
        entry['snc_end'] = toc.snc_time + datetime.timedelta(seconds=(entry['end_sample']-toc.snc_sample)/toc.sample_rate)
        entry['entry_info'] = (self.study_number, self.erd_number, toc.number, -1)
        return entry

    def read_entries(self, entries, t2d=True, buffer=None, interrupt=multiprocessing.Event(), limit=3600):
        if not (self.complete and self.etc.valid):
            return []
        if buffer is None:
            buffer = []
            pro = 4
        else:
            pro =7
        if isinstance(entries, int):
            entries = [entries]
        if len(entries) >= 75:
            with multiprocessing.Pool(processes=pro) as pool:
                results = []
                for num in entries:
                    te = self.etc.etc_data[num]
                    results.append(pool.apply_async(ErdTocFrame, (self.header, self._cf_arr, te, num)))
                for res in results:
                    while (len(buffer) >= limit) and not interrupt.is_set():
                        pass
                    if interrupt.is_set():
                        print('escaping inner')
                        break
                    if t2d:
                        buffer.append(self.toc2dict(res.get(timeout=60)))
                    else:
                        buffer.append(res.get(timeout=60))
        else:
            for num in entries:
                te = self.etc.etc_data[num]
                while (len(buffer) >= limit) and not interrupt.is_set():
                    pass
                if t2d:
                    buffer.append(self.toc2dict(ErdTocFrame(self.header, self._cf_arr, te, num)))
                else:
                    buffer.append(ErdTocFrame(self.header, self._cf_arr, te, num))
        if interrupt.is_set():
            print('got out inner')
        if isinstance(buffer, list):
            return buffer
        else:
            return []

    def read_data(self, etc_file):
        for te in etc_file:
            self._erd_data.append(ErdTocFrame(self.header,
                                              self._cf_arr,
                                              te.offset,
                                              te.samplestamp,
                                              te.sample_span))

        return self._erd_data

    @property  # read once
    def erd_data(self):
        if self._erd_data == []:
            try:
                self.read_data(self.etc)
            except Exception as e:
                print('encountered exception: {:s}'.format(e))

        return self._erd_data

class ErdHeader(object):
    _schema_arrsz = {
        5: 32,
        6: 128,
        7: 1024,
        8: 1024,
        9: 1024,
    }

    def __init__(self, fobj, study_start, study_offset, generic_hdr=None):
        self._fobj = fobj
        if not generic_hdr:
            self._generic_hdr = XLTEK_GenericHeader(self._fobj)
        else:
            self._generic_hdr = generic_hdr

        self.study_start = study_start
        self.study_offset = study_offset
        pos = self._fobj.tell()
        try:
            self.parse()
            self.complete = True
        except struct.error as e:
            self.complete = False
        self._fobj.seek(pos)
        return

    def __getstate__(self):
        state = self.__dict__.copy()
        fobj = state['_fobj']
        state['_fobj'] = (fobj.name, fobj.mode)
        return state

    def __setstate__(self, state):
        name, mode = state['_fobj']
        state['_fobj'] = open(name, mode)
        state['_fobj'].close()
        self.__dict__.update(state)

    def parse(self):
        self._fobj.seek(self._generic_hdr.HDRSIZE)
        self.sample_freq, = struct.unpack('d', self._fobj.read(8))
        self.num_channels, = struct.unpack('I', self._fobj.read(4))
        self.deltabits, = struct.unpack('I', self._fobj.read(4))

        self._arrsz = ErdHeader._schema_arrsz[self._generic_hdr.file_schema]

        fmt = '%dI' % self._arrsz
        var_sz = struct.calcsize(fmt)
        self.phys_chan = list(struct.unpack(fmt, self._fobj.read(var_sz)))

        self.headbox_type = list(struct.unpack('4I', self._fobj.read(16)))
        self.headbox_sn = list(struct.unpack('4I', self._fobj.read(16)))
        self.headbox_sw_version = []
        for hbi in range(4):
            self.headbox_sw_version.append(self._fobj.read(10).decode('ascii').strip('\x00'))

        self.dsp_hw_version = self._fobj.read(10)
        self.dsp_sw_version = self._fobj.read(10)

        self.discardbits, = struct.unpack('I', self._fobj.read(4))

        if self._generic_hdr.file_schema > 7:
            fmt = '%dh' % self._arrsz
            var_sz = struct.calcsize(fmt)
            self.shorted_chan = list(struct.unpack(fmt, self._fobj.read(var_sz)))
            self.freq_factor = list(struct.unpack(fmt, self._fobj.read(var_sz)))

        self._HDRSIZE = self._data_offset = self._fobj.tell()
        return

class ErdTocFrame(object):
    def __init__(self, erd_hdr, cf_arr, toc, num):
        fobj = erd_hdr._fobj
        self.study_start = erd_hdr.study_start
        self.study_offset = erd_hdr.study_offset
        self.sample_rate = erd_hdr.sample_freq

        self.number = num
        self._nbits = erd_hdr.deltabits  # currently only 8 is supported
        self._cf_arr = cf_arr

        self.num_channels = erd_hdr.num_channels
        self._schema_flag = (erd_hdr._generic_hdr.file_schema > 7)

        self.toc_offset = toc.offset
        self.toc_nsamples = toc.sample_span
        self.toc_samplestamp = toc.samplestamp
        self.snc_sample = toc.snc_sample
        self.snc_time = toc.snc_time
        self.parse(fobj)
        return

    def parse(self, fobj):
        self.data = np.zeros((self.toc_nsamples, self.num_channels), dtype='int32')

        # close the fobj if it's open
        tmp_pos = None
        if not fobj.closed:
            tmp_pos = fobj.tell()
            fobj.close()

        # try reading into the data frame using pyxlrd
        try:
           libpyxlrd.read_xlrd_tocframe(self.data, fobj.name, self.toc_offset, self.toc_samplestamp, self.toc_nsamples)
        except Exception as e:
            print('raising')
            raise (e)

        if ((tmp_pos != None) and (fobj.closed)):
            fobj = open(fobj.name)
            fobj.seek(tmp_pos)

        return

class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)

class NormalPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


# Functions #
def filetime_to_datetime(ft):
    # ASSUME LOCAL TIME IS 'America/Los_Angeles'
    timestamp = ft[1]
    timestamp <<= 32
    timestamp |= ft[0]

    utc_datetime = _FILETIME_null_date + datetime.timedelta(microseconds=timestamp / 10)

    return utc_datetime.astimezone(ucsf_time)

def wait(synctime, event=None, interrupt=multiprocessing.Event()):
    interrupt.wait(timeout=synctime)

def file_wait(file, timeout=1):
    if not isinstance(file, pathlib.Path):
        file = pathlib.Path(file)
    flag = True
    then = time.time()
    while not file.is_file():
        if time.time() - then >= timeout:
            flag = False
            break
    return flag

def new_stc_entry(fobj, schema):
    return StcEntry(fobj, schema)

def new_snc_entry(fobj, schema):
    return SncEntry(fobj, schema)

def new_etc_entry(fobj, schema):
    return EtcEntry(fobj, schema)

def make_study(name, direct):
    return XLTEK_Study(name, direct)

def read_erd(erd, entries, t2d=True, buffer=None):
    data = erd.read_entries(entries, t2d, buffer)
    return data

def make_pool(target=None, args=(), info=[]):
    argsq = multiprocessing.Queue()
    returnsq = multiprocessing.Queue()
    args = (argsq,returnsq) + args
    returns = []
    processors =[]
    for data in info:
        argsq.put(data)
    for i in range(0,6):
        processors.append(multiprocessing.Process(target=target,args=args))
        processors[i].start()
    for data in info:
        returns.append(returnsq.get())
    for p in processors:
        p.join()
    return returns

def _make_study_task(argsq, returnsq):
    continue_task = True
    while continue_task:
        if argsq.empty():
            break
        info = argsq.get()
        study = XLTEK_Study(info[0],info[1])
        study.get_etcs()
        study.get_erds()
        print('done')
        returnsq.put(study)

def _buffer_filler(study):
    print('Filling Buffer')
    begin = time.time()
    study.load_buffer()
    print(time.time()-begin)
    print('Done with first fill')
    while not study.stop_events['fill'].is_set():
        study.update_studies()
        study.update_buffer()

########### Main Program ###########
if __name__ == "__main__":
    #ID = 'f69ec446-bb7c-4812-b10a-975b5be7d8bb'
    first = datetime.datetime.now()
    new = XLTEK_Reader('EC164')
    # print(new._studies[-1].sample_rate)
    for study in new._studies:
        print(str(study.study_dir))
        # etc_start = study.stc.stc_data[0].start_stamp - study.sample_offset
        # etc_end = study.stc.stc_data[-1].start_stamp + study.stc.stc_data[-1].sample_span - study.sample_offset
        # start = datetime.timedelta(seconds=(etc_start / study.sample_rate))
        # end = datetime.timedelta(seconds=(etc_end / study.sample_rate))
        # print(study.creation_time+start)
        # print(study.creation_time+end)
        # print(study.sample_offset / study.sample_rate)
    lan = datetime.datetime(2017, 12, 4, 9)
    ent = new.find_entry(lan)
    print(ent)


    #print(datetime.datetime.now()-first)
    # second = datetime.datetime.now()
    # buffer = new.get_to_last(ent)
    # #buffer = new.get_last_data(1)
    # print(datetime.datetime.now() - second)
    # print(len(buffer))
    #second = datetime.datetime.now()
    #buffer = new.get_last_data(3600)
    #print(datetime.datetime.now() - second)
    #second = datetime.datetime.now()
    #buffer = new.update_study_dirs()
    #print(datetime.datetime.now() - second)
    #print(buffer)



    new.start_buffer()
    #buffer, stopper =  new.start_forward(ent)
    #print(buffer)
    #new._buffer_filler()
    print('waiting')
    #time.sleep(10)
    # print(new.subj)
    # print(len(new._buffer))
    study = new._studies[-1]
    while True:
        time.sleep(0.5)
        if len(new._buffer) > 0:
             print(len(new._buffer))
            # print(new._for_buffer[-1].toc_samplestamp)
            #last = new._for_buffer[-1]
            #print(last)
