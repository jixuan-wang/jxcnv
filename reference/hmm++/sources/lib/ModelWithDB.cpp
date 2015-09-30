#include "ModelWithDB.hpp"
#include "Data.hpp"
#include "Exception.hpp"
#include "ModelParams.hpp"
#include "InferredDataFit.hpp"

#include <sstream>
using namespace std;

HMM_PP::ModelWithDB::ModelWithDB(ModelParams* modelParams, string dbFile)
: Model(modelParams), _dbFile(dbFile) {
	// Attach pre-existing dbfile, if exists
	this->attachData();
}

HMM_PP::ModelWithDB::~ModelWithDB() {
	dbShutdown();
}

void HMM_PP::ModelWithDB::attachData() {
	if (!fileExists(_dbFile)) {
		cerr << "Database " << _dbFile << " does not yet exist." << endl;
		return;
	}
	dbCreate();
}

void HMM_PP::ModelWithDB::loadAllData() {
	dbCreate();
	dbClearData();
	HMM_PP::Model::loadAllData();
}

uint HMM_PP::ModelWithDB::addData(Data* d) {
	uint nextIndex = getNumIndividuals();
	dbStoreData(nextIndex, d);
	delete d;

	return nextIndex;
}

const HMM_PP::Data* HMM_PP::ModelWithDB::getData(const uint i) {
	Data* seq = NULL;
	dbRestoreData(i, &seq);
	return seq;
}

HMM_PP::InferredDataFit* HMM_PP::ModelWithDB::getDataFitResult(const uint i) {
	InferredDataFit* fit = NULL;
	dbRestoreDataAndFit(i, &fit);
	return fit;
}

void HMM_PP::ModelWithDB::setDataFitResult(InferredDataFit* idf, const uint ind) {
	cerr << "Storing derived data..." << endl;
	dbStoreFit(ind, idf);
}

void HMM_PP::ModelWithDB::finishData(Data* d) const {
	delete d;
}

void HMM_PP::ModelWithDB::finishDataAndFit(InferredDataFit* idf) const {
	delete idf->getData();
	delete idf;
}

real HMM_PP::ModelWithDB::fit(bool calcDigammas, bool viterbi, bool timeEvents) {
	dbClearFit();
	return HMM_PP::Model::fit(calcDigammas, viterbi, timeEvents);
}

/*
  DB code:
 */
void HMM_PP::ModelWithDB::dbCreate() {
	_datadb.open(_dbFile);

	_datadb.synchronous(false);

	_datadb.query("CREATE TABLE IF NOT EXISTS dat("
			"  indiv_id  INTEGER NOT NULL, "
			"  obs_id    INTEGER NOT NULL, "
			"  val       REAL)");

	_datadb.query("CREATE TABLE IF NOT EXISTS indiv("
			"  indiv_id     INTEGER NOT NULL, "
			"  indiv_label  VARCHAR(20)); ");

	// the 'end-products' we store are the alpha, beta, gamma, and viterbi for each person:
	_datadb.query("CREATE TABLE IF NOT EXISTS alpha("
			"  indiv_id  INTEGER NOT NULL, "
			"  obs_id    INTEGER NOT NULL, "
			"  state_id  INTEGER NOT NULL, "
			"  log10Prob      REAL)");

	_datadb.query("CREATE TABLE IF NOT EXISTS beta("
			"  indiv_id  INTEGER NOT NULL, "
			"  obs_id    INTEGER NOT NULL, "
			"  state_id  INTEGER NOT NULL, "
			"  log10Prob      REAL)");

	_datadb.query("CREATE TABLE IF NOT EXISTS gamma("
			"  indiv_id  INTEGER NOT NULL, "
			"  obs_id    INTEGER NOT NULL, "
			"  state_id  INTEGER NOT NULL, "
			"  log10Prob      REAL)");

	_datadb.query("CREATE TABLE IF NOT EXISTS digamma("
			"  indiv_id  INTEGER NOT NULL, "
			"  obs_id    INTEGER NOT NULL, "
			"  state1_id INTEGER NOT NULL, "
			"  state2_id INTEGER NOT NULL, "
			"  log10Prob      REAL)");

	_datadb.query("CREATE TABLE IF NOT EXISTS viterbi("
			"  indiv_id  INTEGER NOT NULL, "
			"  obs_id    INTEGER NOT NULL, "
			"  state     INTEGER)");


	_stmt_insert_indiv = _datadb.prepare("INSERT OR REPLACE INTO indiv (indiv_id, indiv_label) values (:indiv_id, :indiv_label) ; ");
	_stmt_fetch_indiv = _datadb.prepare("SELECT indiv_label FROM indiv WHERE indiv_id == :indiv_id ; ") ;

	_stmt_get_indiv_n = _datadb.prepare("SELECT COUNT(*) FROM (SELECT DISTINCT indiv_id FROM dat); ");
	_stmt_get_obs_n = _datadb.prepare("SELECT COUNT(*) FROM dat WHERE indiv_id == 0; ");

	_stmt_insert = _datadb.prepare("INSERT OR REPLACE INTO dat (indiv_id, obs_id, val) values(:indiv_id, :obs_id, :val) ; ");
	_stmt_fetch = _datadb.prepare("SELECT val FROM dat WHERE indiv_id == :indiv_id ORDER BY obs_id; ");

	_stmt_fetch_alpha = _datadb.prepare("SELECT obs_id, state_id, log10Prob FROM alpha WHERE indiv_id == :indiv_id ORDER BY obs_id, state_id ; ");
	_stmt_insert_alpha = _datadb.prepare("INSERT OR REPLACE INTO alpha (indiv_id, obs_id, state_id, log10Prob) values(:indiv_id, :obs_id, :state_id, :log10Prob) ; ");
	_stmt_clear_alpha = _datadb.prepare("DELETE FROM alpha WHERE indiv_id == :indiv_id; ");

	_stmt_fetch_beta = _datadb.prepare("SELECT obs_id, state_id, log10Prob FROM beta WHERE indiv_id == :indiv_id ORDER BY obs_id, state_id ; ");
	_stmt_insert_beta = _datadb.prepare("INSERT OR REPLACE INTO beta (indiv_id, obs_id, state_id, log10Prob) values(:indiv_id, :obs_id, :state_id, :log10Prob) ; ");
	_stmt_clear_beta = _datadb.prepare("DELETE FROM beta WHERE indiv_id == :indiv_id; ");

	_stmt_fetch_gamma = _datadb.prepare("SELECT obs_id, state_id, log10Prob FROM gamma WHERE indiv_id == :indiv_id ORDER BY obs_id, state_id ; ");
	_stmt_insert_gamma = _datadb.prepare("INSERT OR REPLACE INTO gamma (indiv_id, obs_id, state_id, log10Prob) values(:indiv_id, :obs_id, :state_id, :log10Prob) ; ");
	_stmt_clear_gamma = _datadb.prepare("DELETE FROM gamma WHERE indiv_id == :indiv_id; ");

	_stmt_fetch_digamma = _datadb.prepare("SELECT obs_id, state1_id, state2_id, log10Prob FROM digamma WHERE indiv_id == :indiv_id ORDER BY obs_id, state1_id, state2_id ; ");
	_stmt_insert_digamma = _datadb.prepare("INSERT OR REPLACE INTO digamma (indiv_id, obs_id, state1_id, state2_id, log10Prob) "
			"values(:indiv_id, :obs_id, :state1_id, :state2_id, :log10Prob) ; ");
	_stmt_clear_digamma = _datadb.prepare("DELETE FROM digamma WHERE indiv_id == :indiv_id; ");

	_stmt_fetch_viterbi = _datadb.prepare("SELECT state FROM viterbi WHERE indiv_id == :indiv_id ORDER BY obs_id ; ");
	_stmt_insert_viterbi = _datadb.prepare("INSERT OR REPLACE INTO viterbi (indiv_id, obs_id, state) values(:indiv_id, :obs_id, :state) ; ");
	_stmt_clear_viterbi = _datadb.prepare("DELETE FROM viterbi WHERE indiv_id == :indiv_id; ");

	//addIndex();
}

void HMM_PP::ModelWithDB::dbShutdown() {
	_datadb.finalize(_stmt_insert_indiv);
	_datadb.finalize(_stmt_fetch_indiv);
	_datadb.finalize(_stmt_get_indiv_n);
	_datadb.finalize(_stmt_get_obs_n);

	_datadb.finalize(_stmt_insert);
	_datadb.finalize(_stmt_fetch);

	_datadb.finalize(_stmt_insert_alpha);
	_datadb.finalize(_stmt_fetch_alpha);
	_datadb.finalize(_stmt_clear_alpha);

	_datadb.finalize(_stmt_insert_beta);
	_datadb.finalize(_stmt_fetch_beta);
	_datadb.finalize(_stmt_clear_beta);

	_datadb.finalize(_stmt_insert_gamma);
	_datadb.finalize(_stmt_fetch_gamma);
	_datadb.finalize(_stmt_clear_gamma);

	_datadb.finalize(_stmt_insert_digamma);
	_datadb.finalize(_stmt_fetch_digamma);
	_datadb.finalize(_stmt_clear_digamma);

	_datadb.finalize(_stmt_insert_viterbi);
	_datadb.finalize(_stmt_fetch_viterbi);
	_datadb.finalize(_stmt_clear_viterbi);

	_datadb.close();
}

uint HMM_PP::ModelWithDB::fetch_n() {
	_datadb.step(_stmt_get_indiv_n);
	uint c = _datadb.get_int(_stmt_get_indiv_n, 0);
	_datadb.reset(_stmt_get_indiv_n);
	return c;
}

uint HMM_PP::ModelWithDB::fetch_obs_n() {
	_datadb.step(_stmt_get_obs_n);
	uint c = _datadb.get_int(_stmt_get_obs_n, 0);
	_datadb.reset(_stmt_get_obs_n);
	return c;
}

void HMM_PP::ModelWithDB::addIndex() {
	_datadb.query("CREATE INDEX IF NOT EXISTS idx1 ON dat(indiv_id, obs_id) ; ");
	_datadb.query("CREATE INDEX IF NOT EXISTS idx2 ON alpha(indiv_id) ; ");
	_datadb.query("CREATE INDEX IF NOT EXISTS idx3 ON beta(indiv_id) ; ");
	_datadb.query("CREATE INDEX IF NOT EXISTS idx4 ON gamma(indiv_id) ; ");
	_datadb.query("CREATE INDEX IF NOT EXISTS idx5 ON digamma(indiv_id) ; ");
	_datadb.query("CREATE INDEX IF NOT EXISTS idx6 ON viterbi(indiv_id); ");
}

void HMM_PP::ModelWithDB::dropIndex() {
	_datadb.query("DROP INDEX IF EXISTS idx1; ");
	_datadb.query("DROP INDEX IF EXISTS idx2; ");
	_datadb.query("DROP INDEX IF EXISTS idx3; ");
	_datadb.query("DROP INDEX IF EXISTS idx4; ");
	_datadb.query("DROP INDEX IF EXISTS idx5; ");
	_datadb.query("DROP INDEX IF EXISTS idx6; ");
}

void HMM_PP::ModelWithDB::dbClearData() {
	_datadb.query("DELETE FROM dat; ");
	_datadb.query("DELETE FROM indiv; ");
	dbClearFit();
}

void HMM_PP::ModelWithDB::dbClearFit() {
	cerr << "Clearing derived information in db..." << endl;
	_datadb.query("DELETE FROM alpha; ");
	_datadb.query("DELETE FROM beta; ");
	_datadb.query("DELETE FROM gamma; ");
	_datadb.query("DELETE FROM digamma; ");
	_datadb.query("DELETE FROM viterbi; ");
}

void HMM_PP::ModelWithDB::dbStoreData(const uint i, HMM_PP::Data* data) {
	// take data from sequence i
	_datadb.begin();
	//dropIndex();

	_datadb.bind_int(_stmt_insert_indiv, ":indiv_id", i);
	_datadb.bind_text(_stmt_insert_indiv, ":indiv_label", data->getId());
	_datadb.step(_stmt_insert_indiv);
	_datadb.reset(_stmt_insert_indiv);

	_datadb.bind_int(_stmt_insert, ":indiv_id", i);
	for (uint j = 0 ; j < data->getNumObservations(); j++) {
		_datadb.bind_int(_stmt_insert, ":obs_id", j);

		data->bindDatapointToDB(_datadb, _stmt_insert, ":val", j);

		_datadb.step(_stmt_insert);
		_datadb.reset(_stmt_insert);
	}

	//addIndex();
	_datadb.commit();
}

void HMM_PP::ModelWithDB::dbRestoreData(const uint i, HMM_PP::Data** seq) {
	// If the object already exists, we can return
	if (*seq)
		return;

	if (!_datadb.attached())
		throw new Exception("No db attached");

	// Otherwise, create a new 'Data' object
	*seq = this->getNewData();

	addIndex();

	_datadb.begin();

	// get ID from database
	_datadb.bind_int(_stmt_fetch_indiv, ":indiv_id", i);
	_datadb.step(_stmt_fetch_indiv);
	(*seq)->setId(_datadb.get_text(_stmt_fetch_indiv, 0));
	_datadb.reset(_stmt_fetch_indiv);

	// Get data-points
	_datadb.bind_int(_stmt_fetch, ":indiv_id", i);
	while (_datadb.step(_stmt_fetch))
		(*seq)->addDatapointFromDB(_datadb, _stmt_fetch, 0);
	_datadb.reset(_stmt_fetch);

	_datadb.commit();
}

void HMM_PP::ModelWithDB::dbStoreFit(const uint i, InferredDataFit* fit)  {
	if (!_datadb.attached())
		return;

	_datadb.begin();
	//dropIndex();

	_datadb.bind_int(_stmt_insert_alpha, ":indiv_id", i);
	for (uint j = 0 ; j < fit->getData()->getNumObservations(); j++)  {
		_datadb.bind_int(_stmt_insert_alpha, ":obs_id", j);
		for (uint state = 0 ; state < getNumHiddenStates() ; state++) {
			_datadb.bind_int(_stmt_insert_alpha, ":state_id", state);
			_datadb.bind_double(_stmt_insert_alpha, ":log10Prob", fit->getAlpha(j,state).getLog10Value());
			_datadb.step(_stmt_insert_alpha);
			_datadb.reset(_stmt_insert_alpha);
		}
	}

	_datadb.bind_int(_stmt_insert_beta, ":indiv_id", i);
	for (uint j = 0 ; j < fit->getData()->getNumObservations(); j++)  {
		_datadb.bind_int(_stmt_insert_beta, ":obs_id", j);
		for (uint state = 0 ; state < getNumHiddenStates() ; state++) {
			_datadb.bind_int(_stmt_insert_beta, ":state_id", state);
			_datadb.bind_double(_stmt_insert_beta, ":log10Prob", fit->getBeta(j,state).getLog10Value());
			_datadb.step(_stmt_insert_beta);
			_datadb.reset(_stmt_insert_beta);
		}
	}

	_datadb.bind_int(_stmt_insert_gamma, ":indiv_id", i);
	for (uint j = 0 ; j < fit->getData()->getNumObservations(); j++)  {
		_datadb.bind_int(_stmt_insert_gamma, ":obs_id", j);
		for (uint state = 0 ; state < getNumHiddenStates() ; state++) {
			_datadb.bind_int(_stmt_insert_gamma, ":state_id", state);
			_datadb.bind_double(_stmt_insert_gamma, ":log10Prob", fit->getGamma(j,state).getLog10Value());
			_datadb.step(_stmt_insert_gamma);
			_datadb.reset(_stmt_insert_gamma);
		}
	}

	if (fit->hasViterbi()) {
		_datadb.bind_int(_stmt_insert_viterbi, ":indiv_id", i);
		for (uint j = 0 ; j < fit->getData()->getNumObservations(); j++)  {
			_datadb.bind_int(_stmt_insert_viterbi, ":obs_id", j);
			_datadb.bind_int(_stmt_insert_viterbi, ":state", fit->getViterbiPath(j));
			_datadb.step(_stmt_insert_viterbi);
			_datadb.reset(_stmt_insert_viterbi);
		}
	}

	if (fit->hasDigamma()) {
		_datadb.bind_int(_stmt_insert_digamma, ":indiv_id", i);
		// diGamma only defined for n-1 consecutive pairs:
		for (uint j = 0; j < fit->getData()->getNumObservations() - 1; j++)  {
			_datadb.bind_int(_stmt_insert_digamma, ":obs_id", j);
			for (uint state1 = 0 ; state1 < getNumHiddenStates() ; state1++) {
				_datadb.bind_int(_stmt_insert_digamma, ":state1_id", state1);
				for (uint state2 = 0 ; state2 < getNumHiddenStates() ; state2++) {
					_datadb.bind_int(_stmt_insert_digamma, ":state2_id", state2);
					_datadb.bind_double(_stmt_insert_digamma, ":log10Prob", fit->getDigamma(j,state1,state2).getLog10Value());
					_datadb.step(_stmt_insert_digamma);
					_datadb.reset(_stmt_insert_digamma);
				}
			}
		}
	}

	//addIndex();
	_datadb.commit();

	// Now we can get rid of this
	fit->clear();
}

void HMM_PP::ModelWithDB::dbRestoreDataAndFit(const uint i, InferredDataFit** fit) {
	// get from database
	if (!_datadb.attached())
		throw new Exception("No db attached");

	// nothing to do?
	if (*fit)
		return;

	HMM_PP::Data* seq = NULL;
	dbRestoreData(i, &seq);
	const uint numObservations = seq->getNumObservations();

	*fit = new InferredDataFit(seq, this);

	addIndex();

	_datadb.begin();

	// Get alphas
	NamedMatrix<real>* alpha = new NamedMatrix<real>();
	alpha->setDims(numObservations, getNumHiddenStates(), 0);
	_datadb.bind_int(_stmt_fetch_alpha, ":indiv_id", i);
	while (_datadb.step(_stmt_fetch_alpha))  {
		const uint t = _datadb.get_int(_stmt_fetch_alpha, 0);
		const uint i = _datadb.get_int(_stmt_fetch_alpha, 1);
		const double v = _datadb.get_double(_stmt_fetch_alpha, 2);
		(*alpha)(t,i) = real(v, true);
	}
	_datadb.reset(_stmt_fetch_alpha);
	(*fit)->setAlpha(alpha);

	// Get betas
	NamedMatrix<real>* beta = new NamedMatrix<real>();
	beta->setDims(numObservations, getNumHiddenStates(), 0);
	_datadb.bind_int(_stmt_fetch_beta, ":indiv_id", i);
	while (_datadb.step(_stmt_fetch_beta))  {
		const uint t = _datadb.get_int(_stmt_fetch_beta, 0);
		const uint i = _datadb.get_int(_stmt_fetch_beta, 1);
		const double v = _datadb.get_double(_stmt_fetch_beta, 2);
		(*beta)(t,i) = real(v, true);
	}
	_datadb.reset(_stmt_fetch_beta);
	(*fit)->setBeta(beta);

	// Get gammas
	NamedMatrix<real>* gamma = new NamedMatrix<real>();
	gamma->setDims(numObservations, getNumHiddenStates(), 0);
	_datadb.bind_int(_stmt_fetch_gamma, ":indiv_id", i);
	while (_datadb.step(_stmt_fetch_gamma))  {
		const uint t = _datadb.get_int(_stmt_fetch_gamma, 0);
		const uint i = _datadb.get_int(_stmt_fetch_gamma, 1);
		const double v = _datadb.get_double(_stmt_fetch_gamma, 2);
		(*gamma)(t,i) = real(v, true);
	}
	_datadb.reset(_stmt_fetch_gamma);
	(*fit)->setGamma(gamma);

	// Get digamma
	InferredDataFit::digamma_t* digamma = new InferredDataFit::digamma_t(numObservations);
	_datadb.bind_int(_stmt_fetch_digamma, ":indiv_id", i);
	while (_datadb.step(_stmt_fetch_digamma))  {
		const uint t = _datadb.get_int(_stmt_fetch_digamma, 0);
		const uint i = _datadb.get_int(_stmt_fetch_digamma, 1);
		const uint j = _datadb.get_int(_stmt_fetch_digamma, 2);
		const double v = _datadb.get_double(_stmt_fetch_digamma, 3);

		if ((*digamma)[t].empty())
			(*digamma)[t].setDims(getNumHiddenStates(), getNumHiddenStates());

		(*digamma)[t](i,j) = real(v, true);
	}
	_datadb.reset(_stmt_fetch_digamma);
	(*fit)->setDigamma(digamma);

	// Get viterbi
	vector<uint>* viterbi = new vector<uint>();
	_datadb.bind_int(_stmt_fetch_viterbi, ":indiv_id", i);
	while (_datadb.step(_stmt_fetch_viterbi))  {
		viterbi->push_back(_datadb.get_int(_stmt_fetch_viterbi, 0));
	}
	_datadb.reset(_stmt_fetch_viterbi);
	if (viterbi->size() == seq->getNumObservations())
		(*fit)->setViterbiPath(viterbi);
	else
		delete viterbi;

	_datadb.commit();
}
