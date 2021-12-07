# xltek_reader

Uses h5py, numpy, and matplotlib


Example:
import xltek_reader
import XLTEK_converter

patient = xltek_reader.XLTEK_Reader("SubjectID", study_ID="XLTEK_ID")

converter = XLTEK_converter.XLTEK_Converter("SubjectID", reader=patient)

converter.start_convert()

converter.stop_convert()
patient.shutdown()