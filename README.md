# Decentralized Power Capping (DPC)
This is the core implementation of DPC agents. In a cluster of servers, DPC agents run on each server which locally computes server power cap such that (i) the
aggregated throughput of the entire cluster is maximized, (ii) workload priorities are taken into account, and (iii) the power usage of the cluster is capped at a certain threshold. 

For more details please visit our publication:

R. Azimi, M. Badiei, X. Zhan, L. Na and S. Reda, "Fast Decentralized Power Capping for Server Clusters", to appear in IEEE Symposium on High-Performance Computer Architecture, 2017.

# Important notes:
* Run DPC as root on all nodes.
* Don’t forget to change the indexes in topo.txt. Each node must have a unique index starting from 0 to number of nodes minus one. 
* If any of the algorithm’s parameters changed such as utility functions etc, for the optimal answer find the optimal epsilon and mu using the provided Matlab file. Then, update these values in DPC_opt function.


# Build instruction
Copy the source code on all the nodes and compile them using make. DPC needs to be run on all the nodes.
Edit the topo.txt budget.txt and workload.txt as input files on all the nodes. (Only index needs to be changed, the rest is the same for all the nodes.)
Run DPC on each node. 
We assumed workload information are known a prior and we select them from the pool of information on runtime. Check out WL_monitor function. 

# Inputs:
topo.txt instruction  
topo.txt has the initial configuration of DPC. The values provided is an example, modify it according to your setup. 

first arg: is the index of the node, index your nodes from 0 to N-1 where N is total number of nodes.
second arg: in total number of nodes (N).  
third arg: the adjacency matrix which indicates the communication topology between DPC agents, its and N by N matrix, if node i is communicating with node j then put 1 in row i column j. The graog from this adjacency matrix only needs to be complete so even a chaine will work.  
fourth: list of the hostnames of your nodes separated by space. In Linux, the hostname and IPs can be defined in /etc/hosts.   
fifth: the base port, all nodes listen to this number plus index to connect to other nodes.  
sixth: number of constraints.  
seven: budget for each constraint.  
eighth: constraint matrix if constraint i has node j put 1 in row i column j.  

budget.txt instruction  
has the total power budget, if the budget changes write the new budget to one of nodes budget.txt and DPC adjust itself.

workload.txt instruction  
write the index of the application and it's priority to this file. In the original implementation, this information was provided by the slurmd but here we read it from a file for simplicity.

# Outputs:
Log.txt: used to log main events and logs of each node.  
DPC_debug.txt: used to log the details of each iteration.   
powerCap.txt: power cap printed every second. In our original implementation, power controller enforce/actuate the cap instead of printing it every second.   

# Contacts
If you use any of our ideas please cite our paper and here is our contact information:  
reza_azimi at brown dot edu or sherif_reda at brown dot edu

