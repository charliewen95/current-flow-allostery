import sqlite3
from sqlite3 import Error
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

def create_new_db(db_file):
    """ create a database connection to an SQLite database """
    conn=None
    try:
        conn=sqlite3.connect(db_file)
        print(sqlite3.version)
    except Exception as e:
        print(e)
    finally:
        if conn:
            conn.close()
    
#Create a Network table in the database and add frame data to it
def create_connection(db_file):
    """ create a connection to the specified SQLite database file (db_file)
        :param db_file: database file
        :return: connection object (or None if not successful)"""
    conn=None
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Exception as e:
        print(e)
        
    return conn

def create_table(conn, table_creation_sql_statement):
    """ create a table using the table_creation_sql_statement
        :param conn: Connection object
        :param table_creation_sql_statement: an sqlite CREATE TABLE statement
        :return:
    """
    try:
        c = conn.cursor()
        c.execute(table_creation_sql_statement)
    except Exception as e:
        print(e)

##################################################
# For step2
### This code adapted from ###
### https://writeonly.wordpress.com/2009/07/16/simple-read-only-sqlalchemy-sessions/
### used to ensure input SQL query string will effectively not have write access
def abort_ro(*args,**kwargs):
    ''' the terrible consequences for trying 
        to flush to the db '''
    print("No writing allowed, tsk!  We're telling mom!")
    return 

def db_setup(connstring,readOnly,echo):
    engine = create_engine(connstring, echo=echo)
    Session = sessionmaker(
        bind=engine, 
        autoflush=not readOnly, 
        autocommit=not readOnly
    )
    session = Session()
    
    if readOnly:
        session.flush = abort_ro   # now it won't flush!

    return session, engine
### ### ###
 
