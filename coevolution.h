#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <limits>
#include <math.h>
#include "lib/rng.h"

using namespace std;

typedef vector<double_t> doubles_t;
typedef vector<int> ints;
typedef pair<int, double> intdouble_t;
typedef std::vector<intdouble_t> intdoubles_t;
#define LOGMSG(m) if (params.writeLog) clog << m << endl
#define tab "\t"
#define newline "\n"
#define INF 1000000000000000

RNG rng;

class Debug{
	public:
	bool set;
	int src;
	int dst;
	int v;
} debug;

struct Model {
	public:
	doubles_t mu;
	doubles_t alpha;
	doubles_t eta;
	doubles_t beta;
	double_t omega_phi;
	double_t omega_kap;
	Model(){};
};

ostream& operator<< (ostream &out, const Model& model){
	out << "w_phi" << tab << "w_kap" << tab << "dummy" << tab << "dummy" << newline;
	out << model.omega_phi << tab << model.omega_kap << tab << -1 << tab << -1 << newline;
	out << "mu" << tab << "alpha" << tab << "eta" << tab << "beta" << newline;
	for (size_t i=0; i < model.mu.size(); i++) {
		out << model.mu[i] << tab << model.alpha[i] << tab << model.eta[i] << tab << model.beta[i] << newline;
	}
	return out;
}

struct ActivityEvent {
	public:
	double t;
	int sid;
	int uid;
	int eventId;
	ActivityEvent* parent_event;
	vector<ActivityEvent*> children;

	ActivityEvent(double tt=0, int ss=0, int uu=0) : t(tt), sid(ss), uid(uu) {}
	friend ostream& operator<< (ostream &out, const ActivityEvent& e);
};

ActivityEvent NullActivityEvent(1000000000.0, -1, -1);

ostream& operator<< (ostream &out, const ActivityEvent& e){
	out  << e.t << "\t" << e.sid << "\t" << e.uid;
	out << tab;
	if (e.parent_event != NULL)
		out << e.parent_event->eventId;
	else
		out << "-1";
	return out;
}

typedef vector<ActivityEvent*> ActivityEvents;

struct NodeEvent {
	double t;
	int uid;
	NodeEvent(double tt=0, int uu=0) : t(tt),  uid(uu) {}
	friend ostream& operator<< (ostream &out, const NodeEvent& e);
};

NodeEvent NullNodeEvent(1000000000.0, -1);

ostream& operator<< (ostream &out, const NodeEvent& e){
	out  <<  e.t << "\t" << e.uid;
	return out;
}

typedef vector<NodeEvent*> NodeEvents;


struct LinkEvent {
	double t;
	int from;
	int to;
	LinkEvent(double tt=0, int ff=0, int oo=0) : t(tt), from(ff), to(oo) {}
	friend ostream& operator<< (ostream &out, const LinkEvent& e);
};

LinkEvent NullLinkEvent(1000000000.0, -1, -1);

ostream& operator<< (ostream &out, const LinkEvent& e){
	out  << e.t << "\t" << e.from << "\t" << e.to;
	return out;
}

typedef vector<LinkEvent*> LinkEvents;

class Network {
	public:
	vector<ints> adj_in;
	vector<ints> adj_out;
	doubles_t creation_time;
	int size;;

	void setSize(int size);

	void addLink(int src, int dst) {
		adj_in.at(dst).push_back(src);
		adj_out.at(src).push_back(dst);
	}

	bool hasLink(int src, int dst) {
		if (adj_out[src].size() < adj_in[dst].size()) {
			for (int j=0; j < adj_out[src].size(); j++)
				if (adj_out[src].at(j) == dst)
					return true;
		} else {
			for (int j=0; j < adj_in[dst].size(); j++)
				if (adj_in[dst].at(j) == src)
					return true;	
		}
		return false;
	}

	friend ostream& operator<< (ostream &out, const Network & network);
};

void Network::setSize(int size) {
	adj_in.resize(size);
	adj_out.resize(size);
	creation_time.resize(size*size);
	this->size = size;
}

ostream& operator<< (std::ostream &out, const Network& network){
	// out << network.size << tab;
	int a [network.size];
	out << newline;
	for (int i=0; i < network.size; i++) {
		for (int j = 0; j < network.size; j++)
			a[j] = 0;
		for (int j=0; j < network.adj_out[i].size(); j++) {
			a[network.adj_out[i][j]] = 1;
		}
		for (int j=0; j < network.size; j++) {
			out << tab << a[j];
		}
		out << newline;
	}
	return out;
}

class History {
	public:
	vector<ActivityEvents*> activityEvents;
	vector<LinkEvent*> linkEvents;

	void setSize(int size) {
		activityEvents.resize(size*size);
		for (int i=0; i < activityEvents.size(); i++)
			activityEvents.at(i) = new ActivityEvents;
	}

	~History() {
		int size = activityEvents.size();
		for (int i=0; i < size; i++)
			delete(activityEvents[i]);
		size = linkEvents.size();
		for (int i=0; i < size; i++)
			delete(linkEvents[i]);
	}
};

struct Params {
	string outputFileName;
	string cascadeFileName;
	string modelFileName;

	double_t T;
	int node_count;
	bool writeLog;
	bool finSp;
	bool random_params;
	double_t sparsity;
	double_t mu_mean;
};

ostream& operator<< (ostream &out, const Params & params){
	out << "T" << tab << "N" << tab << "sp" << tab << "dummy" << newline;
	out  << params.T << tab << params.node_count << tab << params.sparsity << tab << -1;
	return out;
}

Model model;
Network network;
History history;
Params params;

inline double_t phi(double_t t, double_t t_i) {
	return model.omega_phi * exp(-model.omega_phi*(t-t_i));
}


inline double_t kappa(double_t t, double_t t_i) {
	return model.omega_kap * exp(-model.omega_kap*(t-t_i));
}

double eps=1e-8;
static inline int dcmp(double x) {
	return x<-eps?-1:x>eps;
}

size_t sample(doubles_t& d, double sum) {
	double_t u = rng.uniform(0, sum);    
	size_t i=0;
	for(;dcmp(u-d[i])>0 && i<d.size();i++) u-=d[i];
	return (i>=d.size())?d.size()-1:i;
}

class Process {
	protected:
	double_t prev_intensity_time;
	double_t intensity;

	public:
	size_t src;
	size_t dst;
	double_t nextEventTime;

	Process(size_t, size_t);

	virtual double_t getIntensity(double_t time) = 0;

	virtual void updateIntensity(ActivityEvent* activityEvent) = 0;

	virtual  double getSample(double_t curTime) = 0; 

	double_t getNextEventTime() {
		return nextEventTime;
	}

	void computeNextEventTime(double_t curTime) {
		double_t lambda_star = getIntensity(curTime);
		while (1) {
			nextEventTime = curTime + rng.exponential(lambda_star);
			// if (nextEventTime > params.T)
			//	break;
			if (params.finSp) {
				if (params.mu_mean*nextEventTime > params.sparsity)
					break;
			} else {
				if (nextEventTime > params.T)
					break;
			}
			double_t lambda_s = getIntensity(nextEventTime);	
			double_t u = rng.uniform(0, 1);    
			if (u * lambda_star < lambda_s)
				break;
			else
				lambda_star = lambda_s;
		}
	}

};

Process::Process(size_t src, size_t dst) {
	this->src = src;
	this->dst = dst;
}

class ActivityProcess:public Process {
	public:
	ActivityProcess(size_t, size_t);
	~ActivityProcess();
	double_t getIntensity(double_t time);
	void updateIntensity(ActivityEvent* activityEvent);

	double_t getSample(double_t curTime) {
		if (src == dst) {
			return curTime + rng.exponential(model.eta.at(src));
		} else {
			double_t t = curTime;
			double_t I1 = model.beta.at(src)*phi(t, curTime);
			while (1) {
				t += rng.exponential(I1);
				//if (params.mu_mean*t > params.sparsity)
				//	return t;
				
				if (params.finSp) {
					if (params.mu_mean*t > params.sparsity)
						return t;
				} else {
					if (t > params.T)
						return t;
				}

				double_t I2 = model.beta.at(src)*phi(t, curTime);
				double_t u = rng.uniform(0,1);
				if (u * I1 < I2) 
					return t;
				else
					I1 = I2;
			}
		}
	}

};

ActivityProcess::~ActivityProcess() {}

ActivityProcess::ActivityProcess(size_t src, size_t dst):
	Process(src, dst) {
	if (src == dst) 
		intensity =  model.eta[src];
	else
		intensity = 0;
	prev_intensity_time = 0;
}

void ActivityProcess::updateIntensity(ActivityEvent* activityEvent) {
	getIntensity(activityEvent->t);
	if (activityEvent->uid != dst) {
		intensity += model.beta[activityEvent->sid];
	}
}

double_t ActivityProcess::getIntensity(double_t time) {
	if (src == dst)
		intensity -= model.eta[src];

	intensity *= phi(time, prev_intensity_time)/model.omega_phi;

	if (src == dst)
		intensity += model.eta[src];

	prev_intensity_time = time;

	return intensity;	
}


class LinkProcess:public Process {
	public:
	LinkProcess(size_t, size_t);
	~LinkProcess();
	double_t getIntensity (double_t time);
	void updateIntensity(ActivityEvent* activityEvent);
	double_t getSample(double_t curTime) {return -1;}
};

LinkProcess::LinkProcess(size_t src, size_t dst):
	Process(src, dst) {
	if (src != dst) 
		intensity = model.mu[src];

}

LinkProcess::~LinkProcess() {}

double_t LinkProcess::getIntensity(double_t time) {
	if (src == dst || network.hasLink(src,dst))
		return 0;

	intensity -= model.mu[src];
	if (intensity < 0.00000001) {
		prev_intensity_time = time;
		intensity = model.mu[src];
		return intensity;
	}
	intensity *= kappa(time, prev_intensity_time)/model.omega_kap;

	intensity += model.mu[src];

	prev_intensity_time = time;
	return intensity;	
}

void LinkProcess::updateIntensity(ActivityEvent* activityEvent) {
	getIntensity(activityEvent->t);
	intensity += model.alpha[activityEvent->sid];
}

class UpdatableHeap {
	public:
	vector<Process*> a;
	vector<int> heap2id;
	vector<int> id2heap;

	int size;

	UpdatableHeap() {}

	void setSize(int size) {
		this->size = size;
		a.resize(size);
		heap2id.resize(size);
		id2heap.resize(size);
	}

	void swap (int smallest, int i) {
		Process* tmp;
		int tmp_pos;
		tmp = a[smallest];
		a[smallest] = a[i];
		a[i] = tmp;
		tmp_pos = heap2id[smallest];
		heap2id[smallest] = heap2id[i];
		heap2id[i] = tmp_pos;
		id2heap[heap2id[i]] = i;
		id2heap[heap2id[smallest]] = smallest;
	}

	void heapify_down(int i) {
		if (i >= floor(size/2))
			return;
		int l = 2*i+1;
		int r = 2*i+2;

		int smallest = i;
		if (l < size && a[l]->getNextEventTime() < a[smallest]->getNextEventTime())
			smallest = l;
		if (r < size && a[r]->getNextEventTime() < a[smallest]->getNextEventTime())
			smallest = r;
		if (smallest != i) {
			swap(smallest, i);
			heapify_down(smallest);
		}
	}

	void build_heap() {
		for (int i=floor(size/2)-1; i >= 0; i--)
			heapify_down(i);
	}

	int get_min() {
		int min_index = heap2id[0];
		return min_index;
	}

	void heapify_up(int i) {
		if (i==0)
			return;
		int parent = floor((i-1)/2);
		if (a[i]->getNextEventTime() < a[parent]->getNextEventTime()) {
			swap(parent, i);
			heapify_up(parent);
		}
	}

	void update_key(int real_id) {
		int i = id2heap[real_id];
		heapify_up(i);
		heapify_down(i);
	}
};



class NonUpdatableHeap {
	public:
	vector<ActivityEvent*> a;

	int size;

	NonUpdatableHeap
	() {size = 0;}

	void setSize(int size) {
		this->size = size;
		a.resize(size);
	}

	void swap (int smallest, int i) {
		ActivityEvent* tmp;
		tmp = a[smallest];
		a[smallest] = a[i];
		a[i] = tmp;
	}

	void heapify_down(int i) {
		if (i >= floor(size/2))
			return;
		int l = 2*i+1;
		int r = 2*i+2;

		int smallest = i;
		if (l < size && a[l]->t < a[smallest]->t)
			smallest = l;
		if (r < size && a[r]->t < a[smallest]->t)
			smallest = r;
		if (smallest != i) {
			swap(smallest, i);
			heapify_down(smallest);
		}
	}

	void build_heap() {
		for (int i=floor(size/2)-1; i >= 0; i--)
			heapify_down(i);
	}

	ActivityEvent* get_min() {
		return a[0];
	}

	void heapify_up(int i) {
		if (i==0)
			return;
		int parent = floor((i-1)/2);
		if (a[i]->t < a[parent]->t) {
			swap(parent, i);
			heapify_up(parent);
		}
	}

	void insert(ActivityEvent* activityEvent) {
		size = size + 1;
		if (size > a.size()) {
			a.push_back(activityEvent);
		}
		else {
			a.at(size-1) = activityEvent;
		}
		heapify_up(size-1);
	}

	ActivityEvent* exctract_min() {
		ActivityEvent* activityEvent = a.at(0);
		a.at(0) = a.at(size-1);
		size = size -1;
		heapify_down(0);
		return activityEvent;
	}

	void print() {
		for (int i=0; i < size; i++)
			clog << *a[i] << endl;
	}
};