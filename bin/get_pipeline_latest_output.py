# =============================================================
# Author : Fatemeh Dorri
# Created on : June 2020
# =============================================================
'''
Obtain the latest pipeline run for given the run name (for example: '191104_M03829_0264_000000000-CN82P')
'''

import os
import argparse
import pandas as pd
import psycopg2

# details below are for login into database
DATABASE_LOGIN_DETAILS = {}
for databaseName in ["lims", "rd"]:
	DATABASE_LOGIN_DETAILS[databaseName] = {}
	DATABASE_LOGIN_DETAILS[databaseName]["host"]     = "192.168.168.139"
	DATABASE_LOGIN_DETAILS[databaseName]["database"] = databaseName + "_wh"
	DATABASE_LOGIN_DETAILS[databaseName]["user"]     = "bioinfo"
	DATABASE_LOGIN_DETAILS[databaseName]["password"] = "bioinfo"


def connect(LOGIN_DETAILS):  # use, for example, DATABASE_LOGIN_DETAILS["lims"] to log into LIMS
	""" Connect to the PostgreSQL database server """
	conn = None
	try:
		# connect to the PostgreSQL server
		#print('Connecting to the PostgreSQL database...')
		conn = psycopg2.connect(**LOGIN_DETAILS)
		return conn
	except (Exception, psycopg2.DatabaseError) as error:
		print(error)


def get_latest_run(r):
    print("RUNNING {},{}".format(r['Run#'], r['Raw Data']))
    conn   = connect(DATABASE_LOGIN_DETAILS["rd"])
    cursor = conn.cursor()

    SQL_CODE = """SELECT webint_run.rid, webint_runmeta.id, webint_runmeta.run_id, webint_runmeta.key, webint_runmeta.value
                    FROM webint_runmeta
                        INNER JOIN webint_run
                        ON ( webint_runmeta.run_id = webint_run.id )
                    WHERE (webint_runmeta.key = 'is_latest_run'
                        AND webint_run.rid::text LIKE '%{}%'
                        AND webint_runmeta.value = 'True')""".format \
        (r['Raw Data'])

    cursor.execute(SQL_CODE)
    x = cursor.fetchall()
    if len(x) > 1:
        print("WARNING: More than one is tagged true for {}".format(r['Run#'] + ',' + r['Raw Data']))
        print(x)
    elif len(x) == 0:
        print("WARNING: No results is tagged true for {}".format(r['Run#'] + ',' + r['Raw Data']))

    try:
        return x[-1][0]
    except IndexError:
        return ''


def main(args):
    run_list = pd.read_csv(args.run_list, sep=',')
    run_list['latest_run'] = run_list[['Run#', 'Raw Data']].apply(lambda x: get_latest_run(x), axis=1)
    run_list.to_csv(args.output_file, index=False, sep=',')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get the latest pipeline run version')
    parser.add_argument('-r', '--run_list', required=True, help="Path to a csv file with column 'Raw Data' that shows "
                                                                "the list of run_name (for example: "
                                                                "191104_M03829_0264_000000000-CN82P).")
    parser.add_argument('-o', '--output_file', required=True, help="Path to output file")
    args = parser.parse_args()
    main(args)
