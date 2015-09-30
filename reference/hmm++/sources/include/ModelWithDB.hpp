#ifndef __MODEL_WITH_DB_H__
#define __MODEL_WITH_DB_H__

#include "Model.hpp"
#include "sqlwrap.hpp"

namespace HMM_PP {
	class ModelWithDB : public Model {

	public:
		ModelWithDB(ModelParams* modelParams, string dbFile);
		virtual ~ModelWithDB();

	protected:
		virtual uint addData(Data* d);

	public:
		virtual const Data* getData(const uint i);
		virtual InferredDataFit* getDataFitResult(const uint i);

	protected:
		virtual void setDataFitResult(InferredDataFit* idf, const uint ind);

		virtual void loadAllData();

	public:
		virtual void finishData(Data* d) const;
		virtual void finishDataAndFit(InferredDataFit* idf) const;

		// HMM operations
		virtual real fit(bool calcDigammas, bool viterbi, bool timeEvents = false);

		virtual uint getNumIndividuals() { return fetch_n(); }

	private:
		string _dbFile;

		void attachData();

		// storage
		void dbCreate();

		void dbStoreData(const uint, Data*);
		void dbRestoreData(const uint, Data**);

		void dbStoreFit(const uint, InferredDataFit*);
		void dbRestoreDataAndFit(const uint, InferredDataFit**);

		void dbClearData();
		void dbClearFit();

		void dbShutdown();

		uint fetch_n();
		uint fetch_obs_n();
		void addIndex();
		void dropIndex();

		// observed data and gammas and Viterbi solutions
		// (also tmp store for digamma's)
		SQL _datadb;

		sqlite3_stmt* _stmt_fetch;
		sqlite3_stmt* _stmt_insert;
		sqlite3_stmt* _stmt_get_indiv_n;
		sqlite3_stmt* _stmt_get_obs_n;

		sqlite3_stmt* _stmt_fetch_alpha;
		sqlite3_stmt* _stmt_insert_alpha;
		sqlite3_stmt* _stmt_clear_alpha;

		sqlite3_stmt* _stmt_fetch_beta;
		sqlite3_stmt* _stmt_insert_beta;
		sqlite3_stmt* _stmt_clear_beta;

		sqlite3_stmt* _stmt_fetch_gamma;
		sqlite3_stmt* _stmt_insert_gamma;
		sqlite3_stmt* _stmt_clear_gamma;

		sqlite3_stmt* _stmt_fetch_digamma;
		sqlite3_stmt* _stmt_insert_digamma;
		sqlite3_stmt* _stmt_clear_digamma;

		sqlite3_stmt* _stmt_fetch_viterbi;
		sqlite3_stmt* _stmt_insert_viterbi;
		sqlite3_stmt* _stmt_clear_viterbi;

		sqlite3_stmt* _stmt_fetch_indiv;
		sqlite3_stmt* _stmt_insert_indiv;
	};
}

#endif
