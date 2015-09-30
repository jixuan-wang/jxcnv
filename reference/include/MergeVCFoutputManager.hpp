#ifndef __MERGE_VCF_OUTPUT_MANAGER_H__
#define __MERGE_VCF_OUTPUT_MANAGER_H__

#include <string>
#include <list>
#include <algorithm>
#include <map>
#include <vector>
#include <utility>
using namespace std;

#include "GenotypeOutputManager.hpp"

namespace HMM_PP {
	class istreamLineReader;
}

namespace XHMM {

	class MergeVCFoutputManager {

	public:
		MergeVCFoutputManager(list<string>* VCFfiles, list<string>* VCFlistFiles);
		~MergeVCFoutputManager();

		void mergeVCFs(string outFile);

	private:
		typedef list<string> VCFheader;

		typedef map<string, string> InfoMap;
		typedef list<string> InfoKeys;

		typedef map<string, string> Genotypes;

		class Info {
		public:
			Info(InfoMap* infoMap, InfoKeys* infoKeys) : _infoMap(infoMap), _infoKeys(infoKeys) {}

			~Info() {
				if (_infoMap != NULL)
					delete _infoMap;

				if (_infoKeys != NULL)
					delete _infoKeys;
			}

			InfoMap& infoMap() { return *_infoMap; }
			InfoKeys& infoKeys() { return *_infoKeys; }

			const InfoMap& infoMap() const { return *_infoMap; }
			const InfoKeys& infoKeys() const { return *_infoKeys; }

		private:
			InfoMap* _infoMap;
			InfoKeys* _infoKeys;
		};

		class VCFreader;
		list<VCFreader*>* _readers;
		vector<string>* _samples;

		class Variant {
		public:
			Variant(const string& variantRow, const vector<string>* samples);
			~Variant();

			const string& getChr() const { return (*_vals)[GenotypeOutputManager::CHROM]; }
			int getPos();
			const string& getID() const { return (*_vals)[GenotypeOutputManager::ID]; }
			const string& getRef() const { return (*_vals)[GenotypeOutputManager::REF]; }
			const string& getAlt() const { return (*_vals)[GenotypeOutputManager::ALT]; }
			const string& getQual() const { return (*_vals)[GenotypeOutputManager::QUAL]; }
			const string& getFilter() const { return (*_vals)[GenotypeOutputManager::FILTER]; }
			const Info* getInfo();
			const string& getFormat() const { return (*_vals)[GenotypeOutputManager::FORMAT]; }
			const Genotypes* getGenotypes();

			void setInfo(Info* info) { if (info == NULL) return; if (_info != NULL) delete _info; _info = info;}
			void setGenotypes(Genotypes* genotypes) { if (genotypes == NULL) return; if (_genotypes != NULL) delete _genotypes; _genotypes = genotypes;}

		private:
			vector<string>* _vals;

			const vector<string>* _samples;

			int _pos;
			Info* _info;
			map<string, string>* _genotypes;
		};

		class VCFreader {
		public:
			VCFreader(string VCFfile);
			~VCFreader();

			const VCFheader* getHeader() const { return _header; }
			const vector<string>* getSamples() const { return _samples; }

			bool hasNext() { return !_stream->eof(); }
			Variant* next();

		private:
			HMM_PP::istreamLineReader* _stream;

			VCFheader* _header;
			vector<string>* _samples;

			void readHeaderAndSamples();
		};

		class VCFwriter {
		public:
			VCFwriter(HMM_PP::ostreamWriter* out, const VCFheader* header, vector<string>* samples);
			~VCFwriter();

			VCFwriter& operator<<(Variant& v);

			class AlleleStats {
			public:
				AlleleStats(const string& ac, const string& an);
				~AlleleStats();

				const list<uint>* getAC() const { return _ac; }
				uint getAN() const { return _an; }

				AlleleStats& operator+=(const AlleleStats& add);

			private:
				list<uint>* _ac;
				uint _an;
			};

			class AlleleAndNonAlleleInfo {
			public:
				AlleleAndNonAlleleInfo(AlleleStats* allStats, InfoMap* info) : _allStats(allStats), _info(info) {}
				~AlleleAndNonAlleleInfo() {
					if (_allStats != NULL)
						delete _allStats;
					if (_info != NULL)
						delete _info;
				}

				AlleleStats& allStats() { return *_allStats; }
				InfoMap& info() { return *_info; }

				const AlleleStats& allStats() const { return *_allStats; }
				const InfoMap& info() const { return *_info; }

			private:
				AlleleStats* _allStats;
				InfoMap* _info;
			};

			static AlleleAndNonAlleleInfo* splitInfoOnAlleleInfo(const InfoMap& info);

		private:
			HMM_PP::ostreamWriter* _out;
			vector<string>* _samples;
		};

	};

}

#endif
