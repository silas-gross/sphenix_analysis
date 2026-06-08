#! /bin/python3 

import psycopg2 as psql
import sys
import datetime

def main():
    calo = "both"
    start_date=None
    end_date=None
    step=12
    print("Extracting argument")
    if len(sys.argv) > 1:
        calo=sys.argv[1]
        if len(sys.argv) > 2:
            start_date=sys.argv[2]
            if len(sys.argv) > 3:
                end_date=sys.argv[3]
                if len(sys.argv) > 4:
                    step=sys.argv[4]
    if start_date:
        try:
            datetime.strptime(start_date, "%Y-%m-%d")
        except ValueError:
            print("Error: start_date must be in format YYYY-MM-DD")
            return
    
    if end_date:
        try:
            datetime.strptime(end_date, "%Y-%m-%d")
        except ValueError:
            print("Error: end_date must be in format YYYY-MM-DD")
            return
    #connect to the daq database
    try:
        print("testing db connection")
        connection = psql.connect(host="sphenixdaqdbreplica", database="daq")
        cursor = connection.cursor()
    except psql.Error as e:
        return 
    emcal_data=[]
    hcal_data=[]

    if calo in ["emcal", "both"]:
        print("Collecting EMCAL leakage currents...")
        emcal_data = query_leakage_currents(cursor, "emcal_mpodlog", start_date, end_date, emcal)
    
    if calo in ["hcal", "both"]:
        print("Collecting HCAL leakage currents...")
        ohcal_data = query_leakage_currents(cursor, "hcalmpodlog", start_date, end_date, "OHCAL")
        ihcal_data = query_leakage_currents(cursor, "hcalmpodlog", start_date, end_date, "IHCAL")
    write_data_to_file(emcal_data, ihcal_data, ohcal_data, calo, "Leakage_Currents_per_12.txt")
    cursor.close()
    connection.close()
    return

def query_leakage_currents(cursor, table_name, start_date, end_date, calo_type):
        # Build the WHERE clause
    where_conditions = [
        "abs(vmeas - vset) < 1.5",
        "vmeas != 0",
        "vset != 0"
    ]
    read="readtime"
    if "HCAL" in calo_type:
        if "I" in calo_type:
            where_conditions.append("hcal > 23")
        elif "O" in calo_type:
            where_condtions.append("hcal<24")
        read="time"

    # Add date range to WHERE clause if provided
    if start_date and end_date:
        where_conditions.append(f"{read} >= '{start_date}'::date")
        where_conditions.append(f"{read} < '{end_date}'::date + interval '1 day'")
    elif start_date:
        where_conditions.append(f"{read} >= '{start_date}'::date")
    elif end_date:
        where_conditions.append(f"{read} < '{end_date}'::date + interval '1 day'")
    
    where_clause = " AND ".join(where_conditions)
    
    # Build query to average every 12 hours
    query = f"""
    SELECT 
        date_trunc('hour', {read}) + (EXTRACT(hour FROM {read})::int / 12) * interval '12 hours' as time_bin,
        AVG(imeas) as imeas_avg,
        COUNT(*) as sample_count
    FROM {table_name}
    WHERE {where_clause}
    GROUP BY time_bin, sector
    ORDER BY time_bin DESC
    """
    
    try:
        cursor.execute(query)
        data = cursor.fetchall()
        print(f"Retrieved {len(data)} averaged records from {table_name}")
        return data
    except psql.Error as e:
        print(f"Error querying {table_name}: {e}")
        return []

def write_data_to_file(emcal_data, ihcal_data, ohcal_data, calo, filename):
    with open(filename, "w") as file:

        header="time"
        if calo in ["emcal", "both"]:
            header=header+", emcal"
        if calo in ["hcal", "both"]:
            header=header+", ihcal, ohcal"

        file.write(header)
        for i in range(max(len(emcal_data), len(ohcal_data))):
            line=""
            if calo in ["emcal", "both"]:
                line="{emcal_data[i][0]}, {emcal_data[i][1]}"
                if calo in ["both"]:
                    line+=", {ihcal_data[i][1]}, {ohcal_data[i][1]}"
            if calo in ["hcal"]:
                lin="{ihcal_data[i][0]}, {ihcal_data[i][1]}, {ohcal_data[i][1]}"
            file.write(line)
        file.close()
        return
