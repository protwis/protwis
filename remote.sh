#!bin/bash
stamp=$(date "+%Y%m%d%H%M%S")
sql_file=$stamp"_protwis_dump.sql"
echo ">>Writing $sql_file"
pg_dump -U protwis -h localhost -o protwis > $sql_file
echo ">> gzip file $sql_file"
gzip $sql_file
echo ">> DONE"
if [ $1 = "u" ]; then
    echo ">> Upload to test server"
    new_sql_file=$sql_file".gz"
    scp $new_sql_file protwis@138.68.98.19:/protwis/sites/sql_dump.sql.gz
fi

