#!bin/bash
serverName=$(jupyter server list)
portName=$(echo $serverName | cut -d ':' -f 4,4 | cut -d '/' -f 1,1)
PID=$(fuser $portName/tcp | cut -f 2,2)
kill -9 $PID
echo 'killed process $PID, servers still running:'
jupyter server list
 
