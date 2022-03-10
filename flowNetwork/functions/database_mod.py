import sys

#create an empty database
network_database_directory='cx26_GB_Network_Database'
network_database_name='cx26_GB_Network.db'

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
    
create_new_db(network_database_directory+'/'+network_database_name)

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

