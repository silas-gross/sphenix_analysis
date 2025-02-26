import pyodbc

def get_run_numbers(cursor):
    # query = """
    # SELECT runnumber
    # FROM datasets
    # WHERE filename like 'DST_CALO%'
    # GROUP BY runnumber
    # HAVING SUM(events) >= 1000000 AND runnumber > 46619;
    # """
    query = """
    SELECT runnumber
    FROM datasets
    WHERE filename like 'DST_CALO%'
    GROUP BY runnumber
    HAVING SUM(events) >= 500000;
    """
    cursor.execute(query)
    run_numbers = [row.runnumber for row in cursor.fetchall()]
    return run_numbers

def filter_golden_runs(file_catalog_run_numbers, production_cursor):
    query = """
    SELECT runnumber
    FROM goodruns
    """
    production_cursor.execute(query)
    all_good_runs = {row.runnumber for row in production_cursor.fetchall()}

    runs_not_in_goodruns = set(file_catalog_run_numbers) - all_good_runs
    print(f"Number of runs not in the goodruns table: {len(runs_not_in_goodruns)}")

    query = """
    SELECT runnumber
    FROM goodruns
    WHERE (emcal_auto).runclass = 'GOLDEN'
    """
    production_cursor.execute(query)
    golden_runs = {row.runnumber for row in production_cursor.fetchall()}

    # Return the intersection of file_catalog_run_numbers and golden_runs
    return list(golden_runs.intersection(file_catalog_run_numbers))

def main():
    # Connect to the FileCatalog database
    file_catalog_conn = pyodbc.connect("DSN=FileCatalog;UID=phnxrc;READONLY=True")
    file_catalog_cursor = file_catalog_conn.cursor()

    # Get unique run numbers with at least 1 million total events
    file_catalog_run_numbers = get_run_numbers(file_catalog_cursor)
    print(f"Number of runs found in the File Catalog: {len(file_catalog_run_numbers)}")

    # Close the FileCatalog database connection
    file_catalog_conn.close()

    # Connect to the Production database
    production_conn = pyodbc.connect("DSN=Production_write")
    production_cursor = production_conn.cursor()

    # Filter run numbers based on 'GOLDEN' status
    golden_run_numbers = filter_golden_runs(file_catalog_run_numbers, production_cursor)
    golden_run_numbers.sort()

    # Save run numbers to a text file
    runlistfile = 'emcal_goodruns.txt'
    with open(runlistfile, 'w') as f:
        for run_number in golden_run_numbers:
            f.write(f"{run_number}\n")
    print(f"Number of GOLDEN runs saved to {runlistfile}: {len(golden_run_numbers)}")

    # Close the Production database connection
    production_conn.close()

if __name__ == "__main__":
    main()
