# Coevolution
This repository contains the codes for the paper **"Coevolve: A joint point process model for information diffusion and network co-evolution."** 
Mehrdad Farajtabar, Yichen Wang, Manuel Gomez-Rodriguez, Shuang Li, Hongyuan Zha, and Le Song.
In Advances in Neural Information Processing Systems, pp. 1945-1953. 2015.

## COMPILE
To compile run the follwoings: <br>
> g++ -c lib/rng.cpp <br>
> g++ -o coevolution main.cpp rng.o

Then the excutable file "coeovlution" is ready for use.

## RUN
To run with complete  use the following command:
> ./coevolution -N 100 -T 100 -sp 0.004 -finSp 0 -ofn trace.txt -cfn cas.txt -mfn model.txt -wl 0 -mu 0.0001 -alpha 0.5 -eta 0.5 -beta 0.5 -rnd 0 -w_phi 1 -w_kap 1 2> log.txt



## INPUT
The parameters are: <br>
- N:	Number of nodes
- T:	Time limit of the simulation 
- sp: Sparsity of limit of the simulation
- finSp: Finishing with sparsity limit (finsSp=1) or with time limit (finsSp=0)
- ofn: Name of output file containing the trace of activities
- cfn: Name of cascade file containing the statstics of casaces
- mfn: Name of model file containing the parameters of model and simulation
- wl: If wl=1 then log file is created.
- mu: Model parameter for mean of baseline (exogenous) rate for link ceration (c.f. paper)
- alpha: Model parameter for mean of excitory coefficient (indogenous) for link creation (c.f. paper)
- eta: Model parameter for mean of baseline (exogenous) rate for retweet (c.f. paper)
- beta: Model parameter for mean of excitory coefficient (indogenous) for retweet (c.f. paper)
- rnd: If this is set to 1 then the model parameters are set unformly at random with mean specified as above otherwise they are exactly equal to the value specified
- w_phi: The decaying kernel coefficient for link creation
- w_kap: the decaying kernel coefficient for retweet

## OUTPUT
Depending on the input specificaiton you will get up to 4 output files.
- Ouput File (specified by ofn): It contains detailed traces of (link and retweet) events ordered by time of happening. There will be 4 or 5 numbers in each line specified by the following heading:
	type	time	src		dst		parent
	* type: 0 denotes a retweet event and 1 denotes a link event.
	* time: Time of event
	* src: The source node to be retweeted or linked to
	* dst: The node who establishes the link or retweets
	* parent: Exists only for retweet events. It is -1 for the retweets that orginated exgonouesly (actually a tweet) and is set to the number of the event which this tweet is a reshare(retweet) of that one. <br>
- Cascade File (specified by cfn): It contains the statistics of the cascades. More especially, it contains 3 records of data:
 	* Cascade Type: The i-th number in this row contains the number of cascades of type i (Refer to the paper for a specificaton of cascade types)
	* Caccade Depth: The i-th number in this row contains the number of cascades with depth i
	* Cascade Size: The i-th number in this rwo contains the number of cascades of size i (number of nodes in the cascade)
- Model File (specified by mfn): Contains the parameters of model and simulaiton,
		T	N	sp  w_phi	w_kap	
	as specified above.
	Also, then in  N lines it has mu,alpha,eta,beta per node.
- Log File (written when wl=1 and is log.txt): contains a log file of what happens. It will be helpful for develpment.
	
## QUESTIONS
For any question please contact Mehrdad Farajtabar (mehrdad@gatech.edu)
