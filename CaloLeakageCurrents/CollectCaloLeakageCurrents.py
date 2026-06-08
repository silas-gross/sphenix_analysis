import psycopg2 as psql
import sys
from datetime import datetime, 
def main():
    calo = "both"
    if len(sys.argv) > 1:
        calo=sys.argv[1]

    connection = psql.connect(host="sphenixdaqdbreplica", database="daq")
    cursor = connection.cursor()

