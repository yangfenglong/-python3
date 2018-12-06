#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import argparse
import time
from datetime import datetime
from path_monitor import PathMonitor, EventName
from task_auto.novogene import messages_for_logging, logging_tools, global_using
from database import * # for sql
from task_auto.database.mysql_table_models import SampleListPath #for table in sql
import hashlib # for md5



def init_args():
    """
    参数配置函数
    """
    parser = argparse.ArgumentParser(
        description="xjd of projects monitor app")
    parser.add_argument('-p', '--path',
                        help='path of dir to be monitor. Ex:/TJPROJ4/DATA1/HiseqX', required=True)
    parser.add_argument('-o', '--out', help='result file stored dir', required=False,
                        default='./')
    parser.add_argument('-v', '--verbose',
                        help='print more runtime info to help debugging',
                        required=False, action='store_true')
    parser.add_argument('-m', '--mark',
                        help='the special flag in the samplelist path for err operation mark',
                        required=False, default='err_path')
    args = parser.parse_args()
    watched_path = args.path
    out_path = os.path.abspath(args.out)
    is_verbose = args.verbose
    special_path_mark = args.mark

    if is_verbose:
        _LOGGER.info('is_verbose ON')
    else:
        _LOGGER.info('is_verbose OFF')
    return watched_path, out_path, is_verbose,special_path_mark

def get_gongxu_num():
    return 2

def main():
    watched_path, out_path, is_verbose, special_path_mark= init_args()
    events_type = [EventName.GET_IN_CLOSE_WRITE.value,
                   EventName.GET_IN_CREATE.value]

    monitor = PathMonitor(watched_path, events_type, is_verbose)
    if not monitor:
        err_msg = 'Can not create object of PathMonitor. Please check the reason. Exit now.'
        _LOGGER.error(err_msg)
        sys.exit(1)
    monitor.start()
    try:
        #conect sql
        mysql = MySql()
        session=mysql.session
        query = session.query(SampleListPath)
        while True:
            for event_type, list_path in monitor.get_event_watched_path(time_out=1):
                list_path = os.path.abspath(list_path)
                if os.path.isdir(list_path):
                    continue
                if query.filter_by(list_path=list_path).first():
                    file_md5=get_md5(list_path) #funtion in database.py
                    if file_md5 == query.filter_by(list_path=list_path).first().md5:
                        continue
                    else:
                 #       refresh_time=datetime.now()
                        query.filter_by(list_path=list_path).update({"md5": file_md5})
                #        query.filter_by(list_path=list_path).update({"refresh_time": refresh_time})
                else:
                    # id=create_datetime_as_id()
                    newEntry={}
                    newEntry['list_path']=list_path
                    newEntry['create_time']=datetime.now()
                    newEntry['md5']=get_md5(list_path)
                    newEntry['bcl_complete_evalue']='function.not_complete'
                    newEntry['used_flag']='new'
                    newEntry['gongxu']=get_gongxu_num()
                    #if special_path_mark in list_path:
                    #    newEntry['data_force_clean_flag']="yes"
                    #else:
                    #    newEntry['data_force_clean_flag']="no"
                    obj = Dict2cls4smplst(newEntry) # class in database.py
                    session.add(obj) 
                    session.commit() #统一提交，创建数据，在此之前数据库是不会有新增数据的
                if is_verbose:
                    _LOGGER.info('In main {}:{}'.format(event_type, list_path))
                time.sleep(2)        
    except KeyboardInterrupt:
        _LOGGER.info(messages_for_logging.NormalInfo.quit_by_usr_press_ctrl_c('main()'))
    finally:
        session.close()
        monitor.shutdown()

if __name__ == '__main__':
    main()
