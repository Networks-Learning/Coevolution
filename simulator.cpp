#include <iostream>
#include "coevolution.h"

using namespace std;

class Simulator {
	public:
	vector<LinkProcess*> linkProcesses;
	vector<ActivityProcess*> activityProcesses;

	void simulate();
	void analyze_cascades();
	void dfs(ActivityEvent*, int, ints&, int&);
	int event_count;
	Simulator();
	intdouble_t getNextEvent();
	void nextEvents2output();
	void intensities2output(double_t);

	UpdatableHeap linkHeap;
	NonUpdatableHeap activityHeap;
	ActivityEvents allActivityEvents;
};

Simulator::Simulator() {
	LOGMSG(network.size);
	linkProcesses.resize(network.size*network.size);
	activityProcesses.resize(network.size*network.size);
	for (size_t src = 0; src < network.size; src ++) {
		for (size_t dst = 0; dst < network.size; dst ++) {
			int index = src * network.size + dst;
			activityProcesses.at(index) = new ActivityProcess(src, dst);
			// activityProcesses.at(index)->computeNextEventTime(0);
			if (src == dst) {
				double_t t = activityProcesses.at(index)->getSample(0);
				ActivityEvent* activityEvent = new ActivityEvent(t, src, dst);
				activityEvent->parent_event = NULL;
				activityHeap.insert(activityEvent);
			}
		}
	}
	LOGMSG("after activity processes");

	linkHeap.setSize(network.size*network.size);
	for (size_t src = 0; src < network.size; src ++) {
		for (size_t dst = 0; dst < network.size; dst ++) {
			int index = src * network.size + dst;
			linkProcesses.at(index) = new LinkProcess(src, dst);
			linkProcesses.at(index)->computeNextEventTime(0);
			linkHeap.a[index] = linkProcesses.at(index);
			linkHeap.heap2id[index] = index;
			linkHeap.id2heap[index] = index;
		}
	}
	LOGMSG("after link processes");
	linkHeap.build_heap();

	//intensities2output(0);
	//nextEvents2output();
}

void Simulator::nextEvents2output() {
	if (!params.writeLog)
		return;
	for (int i=0; i < linkProcesses.size(); i++) {
		if (i % network.size == 0)
			clog << endl;
		if (i == network.size * network.size)
			clog << endl;
		clog << tab << linkProcesses.at(i)->getNextEventTime();
	}
	clog << newline;
	for (int i=0; i < activityHeap.size; i++)
		clog << *activityHeap.a.at(i) << newline;
	clog << newline;
}

void Simulator::intensities2output(double_t t) {
	if (!params.writeLog)
		return;
	for (int i=0; i < linkProcesses.size(); i++) {
		if (i % network.size == 0)
			clog << newline;
		if (i == network.size*network.size)
			clog << newline;
		clog << tab << linkProcesses.at(i)->getIntensity(t);
	}
	clog << newline;
}

void Simulator::simulate() {
	double_t T = params.T;
	string outputFileName = params.outputFileName;
	string modelFileName = params.modelFileName;

	ofstream modelFile;
	ofstream outputFile;
	outputFile.open(outputFileName.c_str());
	modelFile.open(modelFileName.c_str());
	modelFile << "parameters:" << endl;
	modelFile << params << endl << endl;
	modelFile << "model:" << endl;
	modelFile << model << endl;
	double t = 0;
	int iter = 0;
	double_t cout_step = T/20.0;
	int cout_count = 1;

	int n_lnk = 0, n_twt = 0, n_retwt = 0;

	cout << "time" << tab << tab << "n_twt" << tab << "n_retwt" << tab << "n_lnk" << tab << "sparsity" << endl;
	outputFile << "type" << tab << "time" << tab << "src" << tab << "dst" << tab << "parent" << newline;

	while (1) {
		LOGMSG("\n********************   time = " << t << "    ****    iter = " <<  iter << "   ********************");
		
		//intensities2output(t);
		//nextEvents2output();
		if (t > cout_step * cout_count) {
			cout << "" << t << tab << tab << n_twt << tab << n_retwt << tab << n_lnk << tab << double(n_lnk)/network.size/(network.size-1) << endl;
			cout_count++;
		}
		int index = linkHeap.get_min();
		double_t t_next_link = linkProcesses.at(index)->getNextEventTime();
		double_t t_next_activity = activityHeap.get_min()->t;
		//if (t > params.T)
		//	break;
		if (t_next_activity < t_next_link) {
			t = t_next_activity;
			ActivityEvent* activityEvent = activityHeap.exctract_min();
			size_t src = activityEvent->sid;
			size_t dst = activityEvent->uid;
			outputFile << "0" << tab << *activityEvent << endl;
			LOGMSG(newline << "Activity happned : " << * activityEvent);
			// history.activityEvents[src*network.size + dst]->push_back(activityEvent);
			// processes.at(index)->computeNextEventTime2(t);
			for (int j=0; j < network.adj_out[dst].size(); j++) {
				int v = network.adj_out[dst].at(j);
				double next_t = activityProcesses.at(src*network.size+v)->getSample(t);
				LOGMSG("retweet: " << src << tab << v << tab << params.mu_mean << tab << next_t << tab << params.sparsity);
				if ((params.finSp && params.mu_mean*next_t < params.sparsity) ||
					(!params.finSp && next_t < T)) {
					ActivityEvent* newEvent = new ActivityEvent(next_t, src, v);
					newEvent->parent_event = activityEvent;
					activityHeap.insert(newEvent);
				}
/*				if (! network.hasLink(src, v)) {
					processes.at(network.size*network.size + src*network.size+v)->updateIntensity(activityEvent);
					processes.at(network.size*network.size + src*network.size+v)->computeNextEventTime(t);
					heap.update_key(network.size*network.size + src*network.size+v);
				}
*/
			}
			if (activityEvent->parent_event != NULL) {
				activityEvent->parent_event->children.push_back(activityEvent);
			}
			if (activityEvent->parent_event == NULL) {
				double_t next_t = activityProcesses.at(src*network.size + dst)->getSample(t);
				ActivityEvent* replaceEvent = new ActivityEvent(next_t, src, dst);
				replaceEvent->parent_event = NULL;
				activityHeap.insert(replaceEvent);
			}
		    if (src != dst && !network.hasLink(src, dst)) {
                linkProcesses.at(src*network.size+dst)->updateIntensity(activityEvent);
                linkProcesses.at(src*network.size+dst)->computeNextEventTime(t);
                linkHeap.update_key( src*network.size+dst);
            }
            activityEvent->eventId = n_twt+n_retwt+1;
            allActivityEvents.push_back(activityEvent);
			if (activityEvent->parent_event == NULL) 
				n_twt ++;
			else
				n_retwt++;

		} else {
			t = t_next_link;
			size_t src = index / network.size;
			size_t dst = index % network.size;
			LinkEvent* linkEvent = new LinkEvent(t, src, dst);
			LOGMSG(newline << "Link happened : " << * linkEvent);
			// history.linkEvents.push_back(linkEvent);
			outputFile << "2" << tab << *linkEvent << endl;
			network.adj_in.at(dst).push_back(src);
			network.adj_out.at(src).push_back(dst);
			network.creation_time.at(src*network.size+dst) = t;
			linkProcesses.at(index)->computeNextEventTime(t);
			linkHeap.update_key(index);
			n_lnk++;
		}

		if (params.finSp) {
			if (n_lnk > (params.sparsity*network.size*(network.size-1)))
				break;
		} else {
			if (t > T)
				break;
		}

		if (n_twt + n_retwt > 5000000)
			break;

		iter ++;
	}

	modelFile.close();
	outputFile.close();
}

void Simulator::dfs (ActivityEvent* event, int level, ints& levels, int & max_level) {
	// LOGMSG("dfs : " << "level = " << level << tab << *event);
	levels.push_back(level);
	if (level > max_level)
		max_level = level;
	for (int i=0; i < event->children.size(); i++) {
		dfs(event->children.at(i), level+1, levels, max_level);
	}
}

void Simulator::analyze_cascades() {
	ofstream cascadeFile;
	int cas_count[15];
	int others;
	int max_level_count[params.node_count];
	int cas_size_count[params.node_count];
	for (int i=0; i < params.node_count; i++) {
		max_level_count[i] = 0;
		cas_size_count[i] = 0;
	}
	for (int i=0; i < 15; i++)
		cas_count[i] = 0;
	others = 0;
	cascadeFile.open(params.cascadeFileName.c_str());
	int ii2 = allActivityEvents.size()/1.25;
	for (int i=ii2; i < allActivityEvents.size(); i++) {
		ActivityEvent* event = allActivityEvents.at(i);
		if (event->parent_event != NULL)
			continue;
		ints levels;
		int max_level = 0;
		dfs(event, 0, levels, max_level);
		max_level_count[max_level]++;
		cas_size_count[levels.size()]++;
		sort(levels.begin(), levels.end());
		if (levels.size() == 1){  
			cas_count[0]++;
		} else {
			if (levels.size() == 2) {
				cas_count[1] ++;
			} else {
				if (levels.size() == 3) {
					if (levels.at(2) == 1) {
						cas_count[2]++;	
					} else {
						cas_count[3]++;
					}
				} else {
					if (levels.size() == 4) {
						if (levels.at(3) == 1)  {
							cas_count[4]++;
						} else if (levels.at(3) == 2){
							cas_count[5]++;
						} else {
							cas_count[6]++;
						}
					} else {
						if (levels.size() == 5) {
							if (levels.at(4) == 1) {
								cas_count[7] ++;
							} else {
								if (levels.at(4) == 2) {
									if (levels.at(3) == 1) {
										cas_count[8]++;
									} else {
										if (levels.at(2) == 2) {
											cas_count[11]++;
										} else {
											if (event->children.at(0)->children.size() == 2 ||
												event->children.at(1)->children.size() == 2) {
												cas_count[9]++;
											} else {
												cas_count[10]++;
											}
										}
									}
								} else {
									if (levels.at(4) == 3) {
										if (levels.at(3) == 2) {
											cas_count[12]++;
										} else {
											cas_count[13]++;
										}
									}	else {
										cas_count[14] ++;
									}
								}
							}
						} else {
							others ++;
						}
					}
				}
			}
		}
	}
	
	cascadeFile << "cascades type histogram:" << newline;
	int cut_off = 8;
	for (int i=0;i < cut_off; i++) {
		cascadeFile << cas_count[i] << tab;
	}
	for (int i=cut_off; i < 15; i++)
		others += cas_count[i];
	cascadeFile << others << newline << newline;

	cascadeFile << "cascades depth histogram:" << newline;
	for (int i=0; i < params.node_count; i++)
		cascadeFile << max_level_count[i] << tab;
	cascadeFile << newline << newline;

	cascadeFile << "cascades size histogram:" << newline;
	for (int i=1; i < params.node_count; i++)
		cascadeFile << cas_size_count[i] << tab;
	cascadeFile << newline;
	cascadeFile.close();

	// allocated memory clean-up
	for (int i = 0; i < linkProcesses.size(); i++)
		delete(linkProcesses.at(i));
	for (int i = 0; i < allActivityEvents.size(); i++)
		delete(allActivityEvents.at(i));
		
}



