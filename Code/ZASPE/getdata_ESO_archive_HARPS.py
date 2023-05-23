# HARPS data download script
# Also builds file for running ceres

# =============================================================================
# Setup and imports
# =============================================================================

import sys                                         # argv handling 
import os                                          # file system handling 
import glob
import datetime as dt
from astropy.io import fits


from astroquery.eso import Eso
eso = Eso()
# set eso query row limit - setting to -1 doesn't work properly!
eso.ROW_LIMIT = 100000


# =============================================================================
# General variables
# =============================================================================

# Instrument
Instrument = 'HARPS'

# data directory
Data_Dir = '/home/aaron/Desktop/Second_Earth_Target_Selection/Code/ZASPE/'+Instrument+'/'

# pipeline directory
Pipeline_Dir = '/home/aaron/Desktop/Second_Earth_Target_Selection/Code/ZASPE/ceres/'+Instrument.lower()+'/'

# ESO login
User = 'aaaronmaas '
eso.login(User, store_password=True)

# Calibrations program ID
Calib_ID = '60.A-9036'

# program ID - new to old - inlcudes AKA just in case
Program_ID = ['110.23YQ.001', '0110.C-4071', '109.239V.001', '0109.C-0207',
              '108.22A8.001', '0108.C-0127', 
              '106.21ER.001', '0106.C-0199', '105.20GX.001', '0105.C-0020']

# redownload switch
Redownload = False

# reprocess with ceres switch
Reprocess = True

# File number to retrieve
Max_Rows = 1000

# =============================================================================
# Select between set of nights and single target
# =============================================================================

# request mode selection
Mode = input('Type 0 for set of nights, 1 for single target:\n')
Mode = int(Mode)

# set of dates
if Mode == 0:
    # request start date
    Start_Date = input('Type start date in yyyy mm dd format\n')
    # request end date
    End_Date = input('Type end date in yyyy mm dd format\n')
# target
elif Mode == 1:
    Target = input('Type target name\n')
else:
    print('Mode not understood, exiting')
    sys.exit()

## start date
#Start_Date = '2023 01 01'
#
## end date
#End_Date = '2023 01 02'
#
## get today's date and make string
#Today_Date = dt.date.today()
#Today_string = Today_Date.strftime('%Y %m %d')
#
## Use end date or go to today
#Use_End_Date = True
#
## reset to today if not using end date
#if Use_End_Date == False:
#    End_Date = Today_String 
#

# =============================================================================
# Run the query
# =============================================================================

# initialize list of ceres commands
Ceres_Commands = []

# initialize lists of nights
Nights_Downloaded = []
Nights_Reproc = []
Nights_Skipped = []
Nights_No_Data = []

# Set of nights mode
if Mode==0:
    # convert dates into datetime objects
    Query_Date = dt.datetime.strptime(Start_Date,'%Y %m %d')
    Query_End = dt.datetime.strptime(End_Date,'%Y %m %d')
    
    # loop over nights 
    while Query_Date<=Query_End:
        # initialize query check
        do_query = True
        # check if the night already exists
        if os.path.exists(Data_Dir + Query_Date.strftime('%Y%m%d')):
            if Redownload == False:
                print('Data already exists for '+Query_Date.strftime('%Y %m %d')+', skipping')
                do_query = False
                # save night
                Nights_Skipped.append(Query_Date.strftime('%Y %m %d'))
                # save ceres command if needed
                if Reprocess == True:
                    # store ceres command in list
                    Ceres_Commands.append('python harpspipe.py '+Data_Dir+Query_Date.strftime('%Y%m%d')+' -npools 4\n')

            else:
                print('Data already exists for '+Query_Date.strftime('%Y %m %d')+', moving and redownloading')
                # if there is already an "old" folder, warn and skip
                try:
                    os.rename(Data_Dir + Query_Date.strftime('%Y%m%d'), Data_Dir + Query_Date.strftime('%Y%m%d') + '_old')
                    os.mkdir(Data_Dir + Query_Date.strftime('%Y%m%d'))
                    Nights_Reproc.append(Query_Date.strftime('%Y %m %d'))
                except: 
                    print('Old data already exists for '+Query_Date.strftime('%Y %m %d')+', skipping')
                    do_query = False
                    Nights_Skipped.append(Query_Date.strftime('%Y %m %d'))
                    # save ceres command if needed
                    if Reprocess == True:
                        # store ceres command in list
                        Ceres_Commands.append('python ferospipe.py '+Data_Dir+Query_Date.strftime('%Y%m%d')+' -npools 4\n')

        else:
            print('No data exists for '+Query_Date.strftime('%Y %m %d')+', querying archive')
            os.mkdir(Data_Dir + Query_Date.strftime('%Y%m%d'))
        if do_query:
            # query on program IDs
            # need to loop over old ones, can't really semester-cut because sometines the new ones have been late
            # this assumes only one ID is used per night
            # initialize N_query
            N_Query = 0
            # loop over the IDs
            for i in range(len(Program_ID)):
                Program_Query=eso.query_main(column_filters={'instrument':Instrument,'night':Query_Date.strftime('%Y %m %d'),
                                                             'prog_id':Program_ID[i]})
                # check if there's program data
                try:
                    N_Query = len(Program_Query)
                except:
                    N_Query = 0
                # info message
                print(str(N_Query)+' observations found for Program '+Program_ID[i]+' on '+Query_Date.strftime('%Y %m %d'))
                # if we have found data, end the loop
                if N_Query>0:
                    break
            # if there's program data, get all the data for that night
            if N_Query>0:
                # store the night
                Nights_Downloaded.append(Query_Date.strftime('%Y %m %d'))
                # query calibration data
                Calib_Query=eso.query_main(column_filters={'instrument':Instrument,'night':Query_Date.strftime('%Y %m %d'),
                                                             'prog_id':Calib_ID})
                # download program data
                Data_Files = eso.retrieve_data(Program_Query['Dataset ID'], 
                                               destination=Data_Dir+Query_Date.strftime('%Y%m%d'))
                # download calibration data
                Calib_Files = eso.retrieve_data(Calib_Query['Dataset ID'], 
                                                destination=Data_Dir+Query_Date.strftime('%Y%m%d'))
                # remove binned files while we're at it
                Fits_List = glob.glob(Data_Dir+Query_Date.strftime('%Y%m%d')+'*fits')
                for fitsname in Fits_List:         
                    if fits.getheader(fitsname)['ESO DET WIN1 BINX'] != 1:
                        os.remove(fitsname)
                # store ceres command in list
                Ceres_Commands.append('python harpspipe.py '+Data_Dir+Query_Date.strftime('%Y%m%d')+' -npools 4\n')
            # if no data , save the night name
            else:
                Nights_No_Data.append(Query_Date.strftime('%Y %m %d'))
    
        # move to next night
        Query_Date = Query_Date + dt.timedelta(days=1)

# target mode
elif Mode == 1:
    # initialize query check
    do_query = True
    # check if the target directory already exists
    if os.path.exists(Data_Dir + Target):
        if Redownload == False:
            print('Data already exists for '+Target+', skipping existing nights')
            # do_query = False
        else:
            print('Data already exists for '+Target+', moving and redownloading')
            # if there is already an "old" folder, warn and skip
            try:
                os.rename(Data_Dir + Target, Data_Dir + Target + '_old')
                os.mkdir(Data_Dir + Target)
            except: 
                print('Old data already exists for '+Target+', skipping')
                do_query = False
    else:
        print('No data exists for '+Target+', querying archive') 
        os.makedirs(Data_Dir + Target)
    if do_query:
        # query on target
        # try to resolve on simbad
        Target_Query = eso.query_main(column_filters={'instrument':Instrument, 'target':Target, 'resolver': 'simbad'})
        # check if there's target data
        try:
            N_Query = len(Target_Query)
            Resolver = 'simbad'
        except:
            N_Query = 0
        # if there's none, try to resolve on object keyword
        if N_Query == 0:
            Target_Query=eso.query_main(column_filters={'instrument':Instrument, 'target':Target, 'resolver': 'none'})
            # check if there's target data
            try:
                N_Query = len(Target_Query)
                Resolver = 'none'
            except:
                N_Query = 0        
        # if there's data, download
        if N_Query>0:
            # info message
            print(str(N_Query)+' observations found for Target '+Target+', downloading')
            # sort the dataset IDs
            Dataset_Sort = sorted(Target_Query['Dataset ID'])
            # loop over the files found to identify and download the night
            for i in range(N_Query):
                # reinitialize do_query
                do_query = True
                # get date on file name
                Date_Test = Dataset_Sort[i][6:16]
                # convert into datetime objects
                Query_Date = dt.datetime.strptime(Date_Test,'%Y-%m-%d')
                # see if there's data in that night for the target
                Night_Query = eso.query_main(column_filters={'instrument':Instrument, 'target':Target, 'resolver':Resolver,
                                                             'night':Query_Date.strftime('%Y %m %d')})
                # if there's data, download the target and calib files
                try:
                    Good_Date = len(Night_Query)
                except: 
                    Good_Date = 0
                if Good_Date>0:
                    # save the night
                    Nights_Downloaded.append(Query_Date.strftime('%Y %m %d'))
                    # make directory if it doesn't exist
                    if os.path.exists(Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d')):
                        # skip if not redownloading
                        if Redownload == False:
                            print('Data already exists for '+Query_Date.strftime('%Y %m %d')+', skipping')
                            do_query = False
                            Nights_Skipped.append(Query_Date.strftime('%Y %m %d'))
                            # save ceres command if needed
                            if Reprocess == True:
                                # store ceres command in list
                                Ceres_Commands.append('python harpspipe.py '+Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d')+' -npools 4\n')

                        # otherwise move old data
                        else:
                            print('Data already exists for '+Query_Date.strftime('%Y %m %d')+', moving and redownloading')
                            try:
                                os.rename(Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d'), Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d') + '_old')
                                os.mkdir(Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d'))
                                Nights_Reproc.append(Query_Date.strftime('%Y %m %d'))
                            except: 
                                print('Old data already exists for '+Query_Date.strftime('%Y %m %d')+', skipping')
                                do_query = False
                                Nights_Skipped.append(Query_Date.strftime('%Y %m %d'))
                                # save ceres command if needed
                                if Reprocess == True:
                                    # store ceres command in list
                                    Ceres_Commands.append('python harpspipe.py '+Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d')+' -npools 4\n')

                    else:
                        print('No data exists for '+Query_Date.strftime('%Y %m %d')+', downloading')
                        os.mkdir(Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d'))
                    # check if doing query
                    if do_query:
                        # download program data
                        Data_Files = eso.retrieve_data(Night_Query['Dataset ID'], 
                                                       destination=Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d'))
                        # query calibration data
                        Calib_Query=eso.query_main(column_filters={'instrument':Instrument,'night':Query_Date.strftime('%Y %m %d'),
                                                                     'prog_id':Calib_ID})
                        
                        # download calibration data
                        Calib_Files = eso.retrieve_data(Calib_Query['Dataset ID'], 
                                                        destination=Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d'))
                        # remove binned files while we're at it
                        Fits_List = glob.glob(Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d')+'*fits')
                        for fitsname in Fits_List:         
                            if fits.getheader(fitsname)['ESO DET WIN1 BINX'] != 1:
                                os.remove(fitsname)
                        # store ceres command in list
                        Ceres_Commands.append('python harpspipe.py '+Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d')+' -npools 4\n')
                # if there's no data, go to the previous night
                else:
                    # subtract one from the date
                    Query_Date = Query_Date - dt.timedelta(days=1)
                    # save the night
                    Nights_Downloaded.append(Query_Date.strftime('%Y %m %d'))
                    # make directory if it doesn't exist
                    if os.path.exists(Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d')):
                        # skip if not reprocessing
                        if Redownload == False:
                            print('Data already exists for '+Query_Date.strftime('%Y %m %d')+', skipping')
                            do_query = False
                            Nights_Skipped.append(Query_Date.strftime('%Y %m %d'))
                            # save ceres command if needed
                            if Reprocess == True:
                                # store ceres command in list
                                Ceres_Commands.append('python harpspipe.py '+Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d')+' -npools 4\n')

                        # otherwise move old data
                        else:
                            print('Data already exists for '+Query_Date.strftime('%Y %m %d')+', moving and redownloading')
                            try:
                                os.rename(Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d'), Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d') + '_old')
                                os.mkdir(Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d'))
                                Nights_Reproc.append(Query_Date.strftime('%Y %m %d'))
                            except: 
                                print('Old data already exists for '+Query_Date.strftime('%Y %m %d')+', skipping')
                                do_query = False
                                Nights_Skipped.append(Query_Date.strftime('%Y %m %d'))
                                # save ceres command if needed
                                if Reprocess == True:
                                    # store ceres command in list
                                    Ceres_Commands.append('python harpspipe.py '+Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d')+' -npools 4\n')

                    else:
                        print('No data exists for '+Query_Date.strftime('%Y %m %d')+', downloading')
                        os.mkdir(Data_Dir + Target + '/' + Query_Date.strftime('%Y%m%d'))
                    # check if doing the query
                    if do_query:
                        # redo the night query
                        Night_Query = eso.query_main(column_filters={'instrument':Instrument, 'target':Target, 'resolver':Resolver,
                                                                 'night':Query_Date.strftime('%Y %m %d')})
                        # download program data
                        Data_Files = eso.retrieve_data(Night_Query['Dataset ID'], 
                                                       destination=Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d'))
                        # query calibration data
                        Calib_Query=eso.query_main(column_filters={'instrument':Instrument,'night':Query_Date.strftime('%Y %m %d'),
                                                                     'prog_id':Calib_ID})
        
                        # download calibration data
                        Calib_Files = eso.retrieve_data(Calib_Query['Dataset ID'], 
                                                        destination=Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d'))
                        # remove binned files while we're at it
                        Fits_List = glob.glob(Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d')+'*fits')
                        for fitsname in Fits_List:         
                            if fits.getheader(fitsname)['ESO DET WIN1 BINX'] != 1:
                                os.remove(fitsname)
                        # store ceres command in list
                        Ceres_Commands.append('python harpspipe.py '+Data_Dir+Target+'/'+Query_Date.strftime('%Y%m%d')+' -npools 4\n')


# =============================================================================
# Save CERES commands
# =============================================================================

# write ceres commands to file, and make it executable

# sort the commands
# this also removes duplicates 
Ceres_Commands_Sorted = sorted(list(dict.fromkeys(Ceres_Commands)))


# set of nights
if Mode == 0:
    Ceres_Commands_File = Pipeline_Dir+'run_'+Start_Date.replace(' ','')+'-'+End_Date.replace(' ','')+'.sh'
    f = open(Ceres_Commands_File, 'w')
    for i in range(len(Ceres_Commands_Sorted)):
        f.write(Ceres_Commands_Sorted[i])
    f.close()
    os.system("chmod +x "+Ceres_Commands_File)
elif Mode == 1:
    Ceres_Commands_File = Pipeline_Dir+'run_'+Target+'.sh'
    f = open(Ceres_Commands_File, 'w')
    for i in range(len(Ceres_Commands_Sorted)):
        f.write(Ceres_Commands_Sorted[i])
    f.close()
    os.system("chmod +x "+Ceres_Commands_File)

# print info
if Mode == 0:
    print('Nights downloaded: '+', '.join(night for night in Nights_Downloaded))
    print('Nights reprocessed: '+', '.join(night for night in Nights_Reproc))
    print('Nights skipped: '+', '.join(night for night in Nights_Skipped))
    print('Nights without data: '+', '.join(night for night in Nights_No_Data))
elif Mode == 1:
    print('Nights downloaded for target '+Target+': '+', '.join(night for night in Nights_Downloaded))
    print('Nights reprocessed: '+', '.join(night for night in Nights_Reproc))
    print('Nights skipped: '+', '.join(night for night in Nights_Skipped))

