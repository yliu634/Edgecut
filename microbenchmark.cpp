#include "Ludo/ludo_cp_dp.h"
#include "EdgeCut/ludo_separator.h"

typedef uint32_t K;
typedef uint64_t V;

int main(){

  //Test Cuckoo hashing and Ludo locator.
  {
    LudoSeparator<K, V, 32> ls;
    ls.CommonInit();
    ls.CuckooInit();
    ls.FalsePositiveBenchmark(true);
  }


  //{
  //Test the RTT matrix generation.

  //vector<vector<uint16_t>> la;
  //LatencyMatrix(la, true);


  //prepare user list in mobileuserlist.txt and staticuserlist.txt
  //PreUserList<K, V, 32> ul(true);

  //Test the throughput of the different schemes.
  //LudoSeparator<K, V, 32> ls;
  //ls.TestThroughput(/*Seattle*/false, /*Uniform*/false);
  //}



  /* if you want test Ludo separator with Summary.*/
  //{
  //  LudoSeparator<K, V, 32> ls;   //Summary
  //  ls.Summary();
  //  LudoSeparator<K, V, 32> lp(true); //Ludo Separator
  //  lp.Separator();
  //}

  //std::cout << " " << std::endl;

}





/*
  unordered_map<K, V> baseline;
  K k; V v;
  LFSRGen<K> keyGen(0x1234567801234567ULL, 1000000000, 0);
  LFSRGen<V> valueGen(0x1234567887654321ULL, 100000000, 0);

  while (baseline.size() < 1000000) {

    keyGen.gen(&k);
    valueGen.gen(&v);
    baseline.insert(make_pair(k, v));
  }

  //unordered_map<K, V> baseline;
  ifstream infile;
  infile.open("../userdataset.txt", ios::in);
  assert(infile);

  while(!infile.eof()){
    infile >> k; infile >> v;
    baseline.insert(make_pair(k,v));
  }
  infile.close();


  ofstream outfile;
  outfile.open ("../userdataset.txt");
  assert(outfile);
  for(auto &el : baseline)
  outfile << to_string(el.first) << "\t" << to_string(el.second) <<"\n";
  outfile.close();

 */
/*
 *               for (int i = 0; i < lookupCnt; ++i) {
                  InputBase::distribution = exponential;
                  InputBase::bound = nn;
                  //该如何理解InputBase的存在,掉包吗
                  uint seed = Hasher32<string>()(logName);
                  InputBase::setSeed(seed);

                  uint32_t idx = InputBase::rand();
                  zipfianKeys.push_back(keys[idx]);

 */


