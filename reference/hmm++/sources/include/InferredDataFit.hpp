#ifndef __INFERRED_DATA_FIT_H__
#define __INFERRED_DATA_FIT_H__

#include "utils.hpp"
#include "NamedMatrix.hpp"
#include "params.hpp"
#include "PreciseNonNegativeReal.hpp"

#include <string>
#include <vector>
#include <set>
#include <list>
using namespace std;

class SQL;
struct sqlite3_stmt;

namespace HMM_PP {
	class Data;
	class Model;

	// Inference (fit) results for a particular Data sequence:
	class InferredDataFit {

	public:
		InferredDataFit(const Data* data, Model* model, bool timeEvents = false);
		virtual ~InferredDataFit();

		void clear();

	public:
		inline const Model* getModel() const { return _model; }
		inline const Data* getData() const { return _data; }

		/*
		 * High-level HMM queries
		 */
		real calcLikelihood() const;

	protected:
		real calcLikelihood(const NamedMatrix<real>& alpha, const NamedMatrix<real>& beta, const uint t) const;

	public:
		real calcPosteriorNumerator(const uint t1, const vector<uint>& statePath) const;
		real calcPosterior(const uint t1, const vector<uint>& statePath) const;
		real calcPosteriorComplement(const uint t1, const vector<uint>& statePath) const;

		real calcPosteriorNumerator(const uint t1, const uint t2, const uint state) const;
		real calcPosterior(const uint t1, const uint t2, const uint state) const;
		real calcPosteriorComplement(const uint t1, const uint t2, const uint state) const;

		/*
		 * More sophisticated, high-level HMM queries
		 */
	protected:
		class Constraint;
		real calcLikelihoodGivenConstraints(const uint t1, const uint t2, Constraint* constraints) const;

	public:
		real calcLikelihoodGivenEvent(const uint t1, const vector<uint>& statePath) const;
		real calcLikelihoodGivenEvent(const uint t1, const uint t2, const uint state) const;

		real calcPosteriorNumeratorAssignmentsExcludingStates(const uint t1, const uint t2, const set<uint>& excludeStates) const;
		real calcPosteriorAssignmentsExcludingStates(const uint t1, const uint t2, const set<uint>& excludeStates) const;
		real calcPosteriorComplementAssignmentsExcludingStates(const uint t1, const uint t2, const set<uint>& excludeStates) const;

		/*
		 * Basic HMM operations
		 */
		void forwardBackward(bool calcDigammas);
		void viterbi();

	protected:
		void forwardBackward_Viterbi(NamedMatrix<real>::MarginalizeFunc margFunc, NamedMatrix<real>& alpha, NamedMatrix<real>& beta, NamedMatrix<real>& gamma, bool normalizeGamma);
		void calcGamma(const NamedMatrix<real>& alpha, const NamedMatrix<real>& beta, NamedMatrix<real>& gamma, bool normalize);
		void calcDigamma();

		void calcForwardMessage(const uint t, NamedMatrix<real>& alpha, NamedMatrix<real>::MarginalizeFunc margFunc, const Constraint* constraints = NULL) const;
		void calcBackwardMessage(const uint t, NamedMatrix<real>& beta, NamedMatrix<real>::MarginalizeFunc margFunc, const Constraint* constraints = NULL) const;

		void multiplyMessageTimesEmissionProbs(NamedVector<real>& message, const uint t, const Constraint* constraints = NULL) const;

		class Constraint {
		public:
			Constraint() {}
			virtual ~Constraint() {}

			virtual bool stateIsPermitted(const uint state, const uint t) const = 0;
			bool stateIsForbidden(const uint state, const uint t) const { return !stateIsPermitted(state, t); }
		};

		class HomogeneousConstraint : public Constraint {
		public:
			HomogeneousConstraint(set<uint>* excludedStates) : _excludedStates(excludedStates) {}
			virtual ~HomogeneousConstraint() { delete _excludedStates; }

			bool stateIsPermitted(const uint state, const uint t) const { return _excludedStates->find(state) == _excludedStates->end(); }

			void operator+=(const HomogeneousConstraint& add) { _excludedStates->insert(add._excludedStates->begin(), add._excludedStates->end()); }

		private:
			set<uint>* _excludedStates;
		};

		class TimeDependentConstraint : public Constraint {
		public:
			TimeDependentConstraint() : _timeToConstraints(new map<uint, HomogeneousConstraint*>()) {}
			virtual ~TimeDependentConstraint() {
				for (map<uint, HomogeneousConstraint*>::const_iterator it = _timeToConstraints->begin(); it != _timeToConstraints->end(); ++it)
					delete it->second;
				delete _timeToConstraints;
			}

			bool stateIsPermitted(const uint state, const uint t) const {
				map<uint, HomogeneousConstraint*>::const_iterator findIt = _timeToConstraints->find(t);
				if (findIt == _timeToConstraints->end())
					return true; // no exclusions at time t

				const HomogeneousConstraint* tConstraints = findIt->second;
				return tConstraints->stateIsPermitted(state, t);
			}

			void addConstraintsAtTime(const uint t, HomogeneousConstraint* tConstraints) {
				map<uint, HomogeneousConstraint*>::const_iterator findIt = _timeToConstraints->find(t);
				if (findIt == _timeToConstraints->end())
					(*_timeToConstraints)[t] = tConstraints;
				else {
					(*(*_timeToConstraints)[t]) += (*tConstraints);
					delete tConstraints;
				}
			}

		private:
			map<uint, HomogeneousConstraint*>* _timeToConstraints;
		};

		HMM_PP::InferredDataFit::HomogeneousConstraint* createConstraintsExcludedStates(const set<uint>& excludedStates) const;
		HMM_PP::InferredDataFit::HomogeneousConstraint* createConstraintsOnlyPermitStates(const set<uint>& permittedStates) const;
		set<uint>* getAllStatesList() const;

	public:
		real getAlpha(const uint t, const uint i) const;
		real getBeta(const uint t, const uint i) const;
		real getGamma(const uint t, const uint i) const;

		bool hasDigamma() const;
		real getDigamma(const uint t, const uint i, const uint j) const;

		bool hasViterbi() const;
		uint getViterbiPath(const uint t) const;
		vector<HMM_PP::ullintPair> getViterbiSegments(const set<uint>* excludeSegmentTypes = NULL, const set<uint>* forceSegmentBreakInds = NULL) const;
		void displayViterbiFittedSequence(ostream& stream) const;

		typedef vector<NamedMatrix<real> > digamma_t;

		// Non-const accessor functions (when reading alpha/beta/gamma/viterbi from db)
		inline void setAlpha(NamedMatrix<real>* alpha) { _alpha = alpha; }
		inline void setBeta(NamedMatrix<real>* beta) { _beta = beta; }
		inline void setGamma(NamedMatrix<real>* gamma) { _gamma = gamma; }
		inline void setDigamma(digamma_t* digamma) { _digamma = digamma; }
		inline void setViterbiPath(vector<uint>* vpath) { _vpath = vpath; }

		bool isCalcTimes() const { return _runTimes != NULL; }
		void addCalcTime(const string calc, const uint timeInSec) const;

	protected:
		const Data* _data;
		Model* _model;

		// alpha-probs
		NamedMatrix<real>* _alpha;

		// beta-probs
		NamedMatrix<real>* _beta;

		// _gamma-probs
		NamedMatrix<real>* _gamma;

		// di-gammas
		digamma_t* _digamma;

		// Viterbi path
		vector<uint>* _vpath;

		typedef pair<string, uint> TimedEvent;
		// make mutable so that even const function can change this non-essential feature
		mutable list<TimedEvent>* _runTimes;

#define ADD_SEGMENT_CMD \
		if (excludeSegmentTypes == NULL || excludeSegmentTypes->find(dataVec[curStart]) == excludeSegmentTypes->end()) \
		insertSegments(segments, curStart, t-1, forceSegmentBreakInds);

		template<class DataType>
		static vector<HMM_PP::ullintPair> calcSegments(const vector<DataType>& dataVec, const set<DataType>* excludeSegmentTypes = NULL, const set<uint>* forceSegmentBreakInds = NULL) {
			vector<HMM_PP::ullintPair> segments;
			uint vecLen = dataVec.size();

			if (vecLen > 0) {
				uint curStart = 0;

				uint t = 1;
				for (; t < vecLen; ++t) {
					if (dataVec[t] != dataVec[t-1]) {
						ADD_SEGMENT_CMD;
						curStart = t;
					}
				}
				ADD_SEGMENT_CMD;
			}

			return segments;
		}

		static void insertSegments(vector<HMM_PP::ullintPair>& segments, const uint start, const uint stop, const set<uint>* forceSegmentBreakInds);
	};
}

#endif
