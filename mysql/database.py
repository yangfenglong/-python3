#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据库连接与存储包
"""
import sys
from sqlalchemy import Column, String, create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from datetime import datetime
import hashlib # for md5

_BASE = declarative_base()  # 生成orm基类


class DanZhouBao(_BASE):
    """
    sample.list内容表：t_demultiplexed_sample_list_content
    """
    __tablename__ = 'dandan_zhoubao'

    id = Column(Integer, primary_key=True)
    list_path_id = Column(Integer, ForeignKey(SampleListPath.id), nullable=False)
    flowcell_id = Column(String(100), nullable=False)
    bcl_path = Column(String(200), nullable=False)
    raw_path = Column(String(200), nullable=False)
    lane_id = Column(String(10), nullable=False)
    lib_index = Column(String(20), nullable=False)
    department_id = Column(String(20), nullable=False)
    p7_qc_index = Column(String(20), nullable=False)
    p5_qc_index = Column(String(20), nullable=True)
    building_lib_source = Column(String(40), nullable=False)
    reads_num = Column(BigInteger, nullable=True)
    analysis_state = Column(String(40), nullable=False)
    decision_state = Column(String(40), nullable=False)
    qc_state = Column(String(40), nullable=False, server_default='qc.no')
    species_state = Column(String(40), nullable=False, server_default='species.no')
    create_time = Column(DateTime(timezone=True), nullable=False, server_default=func.now())


    def __repr__(self):
        msg_info = 'department: {}\tlib_index: {}'.format(self.department_id, self.lib_index)
        return msg_info

def get_md5(file):
    md5file=open(file,'rb')
    md5=hashlib.md5(md5file.read()).hexdigest()
    md5file.close()
    return md5

def create_datetime_as_id():
    """
    生成日期字符串用做数据库表记录的主键
    """
    fmt = "%Y%m%d%H%M%S%f"
    return datetime.now().strftime(fmt)

class Dict2cls4Dandan(DanZhouBao):
    """
    convert dict to orm for import to mysqldb
    """
    def __init__(self,dct):
            for key,value in dct.items():
                setattr(self, key, value)

class MySql(object):
    def __init__(self):
        self.session = None
        self.__create_session()

    def __create_session(self):
        """
        建立数据库连接
        """
        # 初始化数据库连接:
        SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://root:mima211314@localhost:3308/mysql'
        engine = create_engine(SQLALCHEMY_DATABASE_URI)

        # 创建DBSession类型:
        DBSession = sessionmaker(bind=engine) # 创建与数据库的会话,这里返回的是个class,不是实例
        if DBSession:
            self.session = DBSession() #生成session实例
            print('success: {}'.format(DBSession))
        else:
            print('error: {}'.format('create fail'))

if __name__ == "__main__":
    MySql()

