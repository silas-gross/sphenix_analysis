#! /bin/python3 

import psycopg2 as psql
import sys
import datetime

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
            where_conditions.append("hcal<24")
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
    query = ""
    if "HCAL" in calo_type:
        query=f"""
        SELECT 
            date_trunc('hour', {read}) + (EXTRACT(hour FROM {read})::int / 12) * interval '12 hours' as time_bin,
            AVG(imeas) as imeas_avg,
            COUNT(*) as sample_count
        FROM {table_name}
        WHERE {where_clause}
        GROUP BY time_bin, hcal
        ORDER BY time_bin ASC
        """
    else:
        query=f"""
        SELECT 
            date_trunc('hour', {read}) + (EXTRACT(hour FROM {read})::int / 12) * interval '12 hours' as time_bin,
            AVG(imeas) as imeas_avg,
            COUNT(*) as sample_count
        FROM {table_name}
        WHERE {where_clause}
        GROUP BY time_bin, sector
        ORDER BY time_bin ASC
        """

    try:
        cursor.execute(query)
        data = cursor.fetchall()
        print(f"Retrieved {len(data)} averaged records from {table_name}")
        avg_data=[]
        holdvalue=-1
        for j in range(len(data)-1):
            if j <= holdvalue:
                j+=1
                continue
            val = data[j][1]
            count = 1
            day = data[j][0].strftime("%x")
            hour = data[j][0].strftime("%H")
            for k in range(j+1,len(data)):
                #print(str(j)+":"+str(k))
                day_k = data[k][0].strftime("%x")
                hour_k = data[k][0].strftime("%H")
                if day == day_k and hour == hour_k:
                    val+=data[k][1]*data[k][2]
                    count+=data[k][2]
              #  print(hour_k)
                if hour_k > hour or day_k > day:
                    holdvalue = k 
                    break
            val=val/float(count)
            avg_data.append([data[j][0],val])
        return avg_data
    except psql.Error as e:
        print(f"Error querying {table_name}: {e}")
        return []

def write_data_to_file(emcal_data, ihcal_data, ohcal_data, calo, filename):
    total_data=[] 
    emcaltime=[emcal_data[x][0] for x in range(len(emcal_data))]
    ihcaltime=[ihcal_data[x][0] for x in range(len(ihcal_data))]
    ohcaltime=[ohcal_data[x][0] for x in range(len(ohcal_data))]
    em_new=[]
    ih_new=[]
    oh_new=[]
    time=[]
    ohcalstart=0
    ihcalstart=0
    if calo in "both":
        for i in range(len(emcaltime)):
            ihcal_pass=False
            ohcal_pass=False
            
            for j in range(ihcalstart, len(ihcaltime)):
                if ihcaltime[j] > emcaltime[i]:
                    ihcalstart=j
                    break
                if ihcaltime[j] == emcaltime[i]:
                    ihcalstart=j
                    ihcal_pass=True
                for k in range(ohcalstart, len(ohcaltime)):
                    if ohcaltime[k] > ihcaltime[j] and ohcaltime[k] <= emcaltime[i]:
                        break
                    if ohcaltime[k] > emcaltime[i]:
                        ohcalstart=k
                        break
                    if ohcaltime[k] == ihcaltime[j] and ohcaltime[k] == emcaltime[i]:
                        ohcalstart=k
                        ohcal_pass=True
                        break
                if ihcal_pass==True and ohcal_pass==True:
                    break
            if ihcal_pass==True and ohcal_pass==True:
                time.append(emcaltime[i])
                em_new.append(emcal_data[i][1])
                ih_new.append(ihcal_data[ihcalstart][1])
                oh_new.append(ohcal_data[ohcalstart][1])
                continue
    elif calo in "hcal":
        for j in range(len(ihcaltime)):
            ohcal_pass = False
            for k in range(ohcalstart, len(ohcaltime)):
                if ohcaltime[k] > ihcaltime[j]:
                    ohcalstart=k
                    break
                if ohcaltime[k] == ihcaltime[j]:
                    ohcalstart=k
                    ohcal_pass=True
                    break
            if ohcal_pass==True:
                time.append(ihcaltime[j])
                ih_new.append(ihcal_data[j][1])
                oh_new.append(ohcal_data[ohcalstart][1])
                continue
    else:
        time=emcaltime
        em_new=[emcal_data[x][1] for x in range(len(emcal_data))]

    with open(filename, "w") as file:

        header="time"
        if calo in ["emcal", "both"]:
            header=header+", emcal"
        if calo in ["hcal", "both"]:
            header=header+", ihcal, ohcal"
        header+="\n"
        file.write(header)
        for i in range(len(time)):
            line=""
            if calo in ["emcal", "both"]:
                line="{0}, {1}".format(time[i], em_new[i])
                if calo in ["both"]:
                    line+=", {0}, {1}".format(ih_new[i], oh_new[i])
            if calo in ["hcal"]:
                line="{0}, {1}, {2}".format(time[i], ih_new[i], oh_new[i])
            line+="\n"
            file.write(line)
        file.close()
        return

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
        datetime.datetime.strptime(start_date, "%Y-%m-%d")
    except ValueError:
        print("Error: start_date must be in format YYYY-MM-DD")

if end_date:
    try:
        datetime.datetime.strptime(end_date, "%Y-%m-%d")
    except ValueError:
        print("Error: end_date must be in format YYYY-MM-DD")
#connect to the daq database
connection=None
cursor=None
try:
    print("testing db connection")
    connection = psql.connect(host="sphnxdaqdbreplica", database="daq")
    cursor = connection.cursor()
    print("connected")
    emcal_data=[]
    ihcal_data=[]
    ohcal_data=[]
    if calo in ["emcal", "both"]:
        print("Collecting EMCAL leakage currents...")
        emcal_data = query_leakage_currents(cursor, "emcal_mpodlog", start_date, end_date, "EMCAL")

    if calo in ["hcal", "both"]:
        print("Collecting HCAL leakage currents...")
        ohcal_data = query_leakage_currents(cursor, "hcalmpodlog", start_date, end_date, "OHCAL")
        ihcal_data = query_leakage_currents(cursor, "hcalmpodlog", start_date, end_date, "IHCAL")
    print("writing to data")
    write_data_to_file(emcal_data, ihcal_data, ohcal_data, calo, "Leakage_Currents_per_12.csv")
    cursor.close()
    connection.close()
except psql.Error as e:
    print(e)

