//
// Created by LIU YI on 2021/6/2.
//

#ifndef EDGECUT_LUDO_SEPARATOR_H
#define EDGECUT_LUDO_SEPARATOR_H

#endif

#pragma once

#include "../utils/hash.h"
#include "../common.h"

#include "../Ludo/ludo_cp_dp.h"
#include "../CuckooPresized/cuckoo_cp_dp.h"
#include "latency_st_pl.h"



template<class K, class V, uint IP = 32>
class LudoSeparator: EdgeCutCommon<K, V, 32>{
public:

    /* length of fingerprint,
     * In original scheme: the value field is composed of : FL + IP (8 + 32)
     * In Separator scheme: the field comsists of : flag bit + SFL + IP (1+8+32)
     * */
    static const uint8_t FL = 10;
    static const uint8_t VL = IP + FL;
    static const uint8_t SFL = 6;
    static const uint8_t SVL = IP + SFL + 1;

    /* length of fingerprint in separator
     * The number of this bloom filter is x and undecided.
     * */
    static const uint8_t ll = 2;
    //static const uint8_t lx;

    //number of mobile and stationary users
    using EdgeCutCommon<K, V, 32>::nn;
    using EdgeCutCommon<K, V, 32>::mm;
    //static const K nn = 100000;
    //static const K mm = 5 * nn;

    //value field length = fingerprint length + IP length;
    static const uint64_t FPMask = (1ULL << FL) - 1;
    static const uint64_t ValueMask = (1ULL << VL) - 1;

    static const uint64_t SFPMask = (1ULL << SFL) - 1;
    static const uint64_t SValueMask = (1ULL << SVL) - 1;

    static const uint64_t IPMask = (1ULL << IP) - 1;
    static const uint64_t fallbackMask = (1ULL << ll) - 1;

    ControlPlaneLudo<K, V, VL> cp;
    ControlPlaneLudo<K, V, SVL> scp;
    ControlPlaneCuckoo<K, V, uint32_t, true, 8> rp;

    DataPlaneLudo<K, V, VL> *dp;
    DataPlaneLudo<K, V, SVL> *sdp;
    DataPlaneCuckoo <K, V, uint32_t, 8> *rdp;

    vector<vector<uint16_t>> la;

    vector<K> keys;
    vector<V> values;

    uint16_t ns;


    LudoSeparator(){

    }

    explicit LudoSeparator(const bool flag){
      string str = "../userdataset.txt";
      if (flag)
        CommonInit(str, true);
      else
        CommonInit(str);
    }

    void CommonInit(){
      LFSRGen<K> keyGen(0x1234567801234567ULL, 10 * nn, 0);
      LFSRGen<V> valueGen(0x1234567887654321ULL, nn, 0);

      keys.resize(nn);
      values.resize(nn);
      cp.resizeCapacity(nn, false);
      for (uint32_t i = 0; i < nn; i ++) {
        keyGen.gen(&keys[i]);
        V v;
        valueGen.gen(&v);
        // values[i] = v & ValueMask;
        values[i] = ((FastHasher64<K>(123)(keys[i]) << 32) + (v & IPMask)) & ValueMask;
        cp.insert(keys[i], values[i]);
      }
      dp = new DataPlaneLudo<K, V, VL> (cp);
      //V vt = dp->lookUp(keys[3]);
    }

    // Init Cuckoo Summary
    void CuckooInit(){
      LFSRGen<K> keyGen(0x1234567801234567ULL, 10 * nn, 0);
      LFSRGen<V> valueGen(0x1234567887654321ULL, nn, 0);

      keys.resize(nn);
      values.resize(nn);
      rp.Clear(nn);

      for (uint32_t i = 0; i < nn; i ++) {
        keyGen.gen(&keys[i]);
        V v;
        valueGen.gen(&v);
        // values[i] = v & ValueMask;
        values[i] = ((FastHasher64<K>(123)(keys[i]) << 32) + (v & IPMask)) & ValueMask;
        rp.insert(keys[i], values[i]);
      }
      rdp = new DataPlaneCuckoo <K, V, uint32_t, 8> (rp);

      // V vt; rdp->lookUp(4, vt);
    }

    /* for no separator ludo co and dp*/
    void CommonInit(const string & str){
      cp.resizeCapacity(nn, false);
      //unordered_map<K, V> baseline;
      uint32_t i = 0; K k; V v;
      ifstream infile;
      infile.open(str, ios::in);
      assert(infile);
      while(!infile.eof() && i < nn){
        infile >> k; infile >> v;
        v = ((FastHasher64<K>(123)(k) << 32) + (v & IPMask)) & ValueMask;
        assert(v <= (1ULL << VL));
        cp.insert(k, v);
        i ++;
      }
      infile.close();

      dp = new DataPlaneLudo<K, V, VL> (cp);
    }

    /* for ludo separator with ll length additional fingerprint and dp*/
    void CommonInit(const string & str, const bool flag){

      //unordered_map<K, V> baseline;
      scp.resizeCapacity(nn, false);
      uint32_t i = 0; K k; V v;
      ifstream infile;
      infile.open(str, ios::in);
      assert(infile);
      while(!infile.eof() && i < nn){
        infile >> k; infile >> v;
        v = ((FastHasher64<K>(123)(k) << 32) + (v & IPMask)) & ( (1ULL<< (SVL-1)) - 1);
        assert(v < (1ULL << SVL));
        scp.insert(k, v);
        i ++;
      }
      infile.close();

      sdp = new DataPlaneLudo<K, V, SVL> (scp);
    }

    void FalsePositiveBenchmark(const bool CuckooTable){
      LFSRGen<K> keyGen(0x1234567801234567ULL, 10 * nn, 0);
      LFSRGen<V> valueGen(0x1234567887654321ULL, nn, 0);

      unordered_map<K, V> mp;
      K k; V v, vt;
      uint32_t lcount(0), rcount(0);

      for (uint32_t i = 0; i < nn; i ++) {
        keyGen.gen(&k);
        valueGen.gen(&v);
        mp.insert(make_pair(k,v));
      }

      //Test false positive of rdp and dp;
      for(const auto el: mp){
        k = el.first; v = el.second;
        rdp->lookUp(4, vt);
        if (vt != v)
          rcount ++;
        vt = dp->lookUp(k);
        if (vt != v)
          lcount ++;
      }

      std::cout << to_string(rcount) << " : " << lcount << std::endl;

      // Test two locators' lookUp throughput
      Clocker rpt("Cuckoo lookup");
      for (int i = 0; i < nn; i ++) {
        keyGen.gen(&k);
        rdp->lookUp(k, vt);
      }
      rpt.stop();
      Clocker dpt("Ludo lookup");
      for (int i = 0; i < nn; i ++) {
        keyGen.gen(&k);
        dp->lookUp(k);
      }
      dpt.stop();


      // Test two locators' Update throughput
      Clocker zqq("Cuckoo Update");
      for (auto const &el: mp) {
        rp.changeValue(el.first, 3);
      }
      zqq.stop();

      Clocker ly("Ludo Update");
      for (auto const &el: mp) {
        cp.changeValue(el.first, 3);
      }
      ly.stop();

      // Test Meomory
      uint64_t memcuckoo = rdp->getMemoryCost();
      uint64_t memludo = dp->getMemoryCost();
      std::cout << memcuckoo << " : " << memludo << std:: endl;


    }

    void FalsePositiveBenchmark(){
      uint32_t mm(20000000), count(0);
      LFSRGen<K> keyGen(0x1234588801204569ULL, 10*nn, 0);
      vector<K> testKeys;
      testKeys.resize(mm);
      for (uint32_t i = 0; i < mm; i ++) {
        keyGen.gen(&testKeys[i]);
      }
      for (uint32_t i = 0; i < mm; i ++) {
        V v = dp->lookUp(testKeys[i]);
        assert(v <= 1ULL << VL);
        if ((FastHasher64<K>(123)(testKeys[i]) & FPMask) == (((v & ValueMask) >> 32) & FPMask))
          count ++;
      }

      std::cout << "hello count:\n" << count << endl;
    }

    void FalsePositiveBenchmark(const string &str){
      uint32_t count(0), i(0);
      K k; V v;
      ifstream infile;
      infile.open(str, ios::in);
      assert(infile);
      while(!infile.eof() && i < (nn + mm)){
        i ++; infile >> k; infile >> v;
        if (i < nn)
          continue;
        v = dp->lookUp(k);
        assert(v <= 1ULL << VL);
        if ((FastHasher64<K>(123)(k) & FPMask) == (((v & ValueMask) >> 32) & FPMask))
          count++;
      }
      infile.close();
      std::cout << "Original fingerprint count:" << count << endl;
    }

    void fallbackTable(unordered_map<uint64_t, uint8_t> &fallback, const string &str){
      uint32_t count(0), i(0), bid(0);
      uint8_t sid(0);
      K k; V v;

      ifstream infile;
      infile.open(str, ios::in);
      assert(infile);
      while(!infile.eof() && i < (nn + mm)){
        i ++; infile >> k; infile >> v;
        if (i < nn)
          continue;
        sdp->lookUp(k, v, bid, sid);
        //assert(v < (1ULL << SVL));
        if ((FastHasher64<K>(123)(k) & SFPMask) == (((v & ValueMask) >> 32) & SFPMask)){
          count ++;
          sdp->applyUpdate(bid, sid, v + (1ULL << (SVL-1)));
          fallback.insert(make_pair( ((bid << 2) + sid)&((1ULL << 34)-1), FastHasher64<K>(1)(k) & fallbackMask));
        }

      }
      infile.close();
      std::cout << "fallback size count: "<< fallback.size()<<"/"<<count << endl;

    }

    void FalsePositiveBenchmark(const unordered_map<uint64_t, uint8_t> & fallback, const string & str){
      uint32_t count(0), i(0), bid(0);
      std::unordered_map<uint64_t, uint8_t>::const_iterator got;
      uint8_t sid(0);
      K k; V v;

      ifstream infile;
      infile.open(str, ios::in);
      assert(infile);
      while(!infile.eof() && i < (nn + mm)){
        i ++; infile >> k; infile >> v;
        sdp->lookUp(k, v, bid, sid);
        if(i < nn){
          if (v >= 1ULL << SVL){
            got = fallback.find(((bid << 2) + sid) & ((1ULL << 34)-1) );
            assert(got != fallback.end());
            if(got->second != (FastHasher64<K>(1)(k) & fallbackMask))
              count ++;
            }
          } else {
          if (v >= 1ULL << SVL){
            got = fallback.find(((bid << 2) + sid) & ((1ULL << 34)-1) );
            assert(got != fallback.end());
            if(got->second == (FastHasher64<K>(1)(k) & fallbackMask))
              count ++;
          }
        }
      }
      infile.close();
      std::cout << "Ludo Separator fp count:" << count << std::endl;
   }

    //**************************************************

    void Insert_cp_dp(const string &str){
      CommonInit(str);
    }

    void Insert_scp_sdp(const string &str){
      CommonInit(str, true);
    }

    void fallbackTableGen(unordered_map<uint64_t, uint8_t> &fallback, const string &str){
      fallback.clear();
      uint32_t count(0), i(0), bid(0);
      uint8_t sid(0);
      K k; V v;

      ifstream infile;
      infile.open(str, ios::in);
      assert(infile);
      while(!infile.eof() && i < mm){
        i ++; infile >> k; infile >> v;
        sdp->lookUp(k, v, bid, sid);
        //assert(v < (1ULL << SVL));
        if ((FastHasher64<K>(123)(k) & SFPMask) == (((v & ValueMask) >> 32) & SFPMask)){
          count ++;
          sdp->applyUpdate(bid, sid, v + (1ULL << (SVL-1)));
          fallback.insert(make_pair( ((bid << 2) + sid)&((1ULL << 34)-1), FastHasher64<K>(1)(k) & fallbackMask));
        }

      }
      infile.close();
      std::cout << "\nfallback size count: "<< fallback.size() << "/" << count << "\n" << endl;

    }

    void TestThroughput(const bool Seattle = true, const bool Uniform = true) {
      LatencyMatrix(la, Seattle);
      ns = la.size();

      vector<K> TestKeys;
      vector<V> TestValues;
      uint32_t lookUpCnt = 10000;

      bool StaticUsers = true;

      unordered_map<uint64_t, uint8_t> fallback;

      /* ludo hashing */
      Insert_cp_dp("../mobileuserlist.txt");
      Insert_scp_sdp("../mobileuserlist.txt");

      // fallback separator table generation and train
      // using first m static users in staticuserlist.txt file.
      fallbackTableGen(fallback, "../staticuserlist.txt");

      // Test keys generation
      // We can generate 2 kinds of Test users keys
      // 1). only mobile users will be generated. (4 args)
      // 2). N mobile users + \alphe * N static users.(5 args),
      // p.s. nnmm is the separating point of mobile users & static users.
      PreUserList<K, V, 32> pl(Seattle);
      pl.TestKeysGen1(TestKeys, TestValues, lookUpCnt, Uniform);
      // pl.TestKeysGen(TestKeys, TestValues, lookUpCnt, uniform, /*static users contained */ true, nnmm);


      /* Throughput Test (no static users contained) */

      vector<uint32_t> RandomNoSeparator, ShortCutNoSeparator, CloudServerOnly, ProactiveReplicas;

      for(uint32_t i = 1; i < ns; i ++){

        for(uint32_t j = 0; j < lookUpCnt ; j ++) {

          V rl = dp->lookUp(TestKeys[j]);
          rl = rl & IPMask;
          srand(time(NULL));

          if ((rand() % 2) == 0) {
            RandomNoSeparator.push_back(la[i][rl]);
          } else {
            RandomNoSeparator.push_back(la[i][0]);
          }

          if(la[i][rl] <= la[i][0]){
            ShortCutNoSeparator.push_back(la[i][rl]);
          } else{
            ShortCutNoSeparator.push_back(la[i][0]);
          }

          CloudServerOnly.push_back(la[i][0]);
          ProactiveReplicas.push_back(la[i][0] + rand() % 10);
        }

      }

      StatisticsForTraffic(RandomNoSeparator, "Ludo + Random + No Separator");
      StatisticsForTraffic(ShortCutNoSeparator, "Ludo + ShortCut + No Separator");
      StatisticsForTraffic(CloudServerOnly, "Cloud Server Only Approach");
      StatisticsForTraffic(ProactiveReplicas, "Proactive 2 Replicas Approach");



      uint32_t bid(0), nnmm(0); uint8_t sid(0); V srl(0);

      nnmm = (uint32_t) (round(lookUpCnt / ((mm/nn) + 1)));


      if(StaticUsers){

        RandomNoSeparator.clear();
        ShortCutNoSeparator.clear();
        CloudServerOnly.clear();
        ProactiveReplicas.clear();

        vector<uint32_t> ShortCutYesSeparator;
        std::unordered_map<uint64_t, uint8_t>::const_iterator got;
        // if you want to generate Testkays, please clear your vectors.
        pl.TestKeysGen2(TestKeys, TestValues, lookUpCnt, uniform, /*static users contained */ true);

        for (uint32_t i = 1; i < ns; i ++){

          // for mobile users part
          for(uint32_t j = 0; j < nnmm; j ++){

            V rl = dp->lookUp(TestKeys[j]);
            sdp->lookUp(TestKeys[j], srl, bid, sid);

            rl = rl & IPMask;
            srand(time(NULL));

            // choose one server randomly.
            if ((rand() % 2) == 0) {
              RandomNoSeparator.push_back(la[i][rl]);
            } else {
              RandomNoSeparator.push_back(la[i][0]);
            }
            // choose near one using EdgeCut.
            if(la[i][rl] <= la[i][0]){
              ShortCutNoSeparator.push_back(la[i][rl]);
            } else{
              ShortCutNoSeparator.push_back(la[i][0]);
            }

            CloudServerOnly.push_back(la[i][0]);
            ProactiveReplicas.push_back(la[i][0] + rand() % 10);

            if(srl < (1ULL << SVL)){
              srl = srl & IPMask;
              ShortCutYesSeparator.push_back(la[i][srl] <= la[i][0]? la[i][srl] : la[i][0]);
            } else {
              srl = srl & IPMask;
              got = fallback.find(((bid << 2) + sid) & ((1ULL << 34)-1) );
              assert(got != fallback.end());
              if(got->second == (FastHasher64<K>(1)(TestKeys[j]) & fallbackMask))
                ShortCutYesSeparator.push_back(la[i][0]);
              else
                ShortCutYesSeparator.push_back(la[i][srl] <= la[i][0]? la[i][srl] : la[i][0]);
            }


          }

          for(uint32_t j = nnmm; j < lookUpCnt; j++){

            V rl = dp->lookUp(TestKeys[j]);
            sdp->lookUp(TestKeys[j], srl, bid, sid);

            //if the static users was justified as mobile users.
            if ((FastHasher64<K>(123)(TestKeys[j]) & FPMask) == (((rl & ValueMask) >> 32) & FPMask)) {
              uint16_t tmp = la[i][0];
              rl = rl & IPMask;
              if (la[i][rl] < la[i][0])
                tmp += la[i][rl];
              RandomNoSeparator.push_back(tmp);
              ShortCutNoSeparator.push_back(tmp);
            } else{
              // rl = rl & IPMask;
              RandomNoSeparator.push_back(la[i][0]);
              ShortCutNoSeparator.push_back(la[i][0]);
            }

            CloudServerOnly.push_back(la[i][0]);
            ProactiveReplicas.push_back(la[i][0] + rand() % 10);

            if(srl < (1ULL << SVL)){
              srl = srl & IPMask;
              ShortCutYesSeparator.push_back(la[i][0]);
            } else {
              srl = srl & IPMask;
              got = fallback.find(((bid << 2) + sid) & ((1ULL << 34)-1) );
              assert(got != fallback.end());
              if(got->second == (FastHasher64<K>(1)(TestKeys[j]) & fallbackMask))
                ShortCutYesSeparator.push_back(la[i][0]);
              else {
                uint16_t tmp = la[i][0];
                if (la[i][srl] <= la[i][0])
                   tmp += la[i][srl];
                ShortCutYesSeparator.push_back(tmp);
              }
            }

          }

        }

        std::cout << "\n" << std::endl;
        StatisticsForTraffic(RandomNoSeparator, "Ludo + Random + No Separator");
        StatisticsForTraffic(ShortCutNoSeparator, "Ludo + ShortCut + No Separator");
        StatisticsForTraffic(ShortCutYesSeparator, "Ludo + ShortCut + Yes Separator");
        StatisticsForTraffic(CloudServerOnly, "Cloud Server Only Approach");
        StatisticsForTraffic(ProactiveReplicas, "Proactive 2 Replicas Approach");


      } //end if(Seattle)


    }//end Test Throughput

    void StatisticsForTraffic(const vector<uint32_t> &list, const string &str){
      double sum(0), CrossTraffic(0);

      for(auto const &el: list){
        sum += el;
        if (el >= 100)
          CrossTraffic ++;
      }
      double mean = sum/list.size();
      double CrossTrafficRate = (CrossTraffic/list.size()) * 100;
      std::cout << str << ": " << mean << " " << round(CrossTrafficRate) << "%" << endl;

    }


    void Summary(){
      string str = "../userdataset.txt";
      CommonInit(str);
      FalsePositiveBenchmark(str);
    }

    void Separator(){
      string str = "../userdataset.txt";
      unordered_map <uint64_t, uint8_t> fallback;
      fallbackTable(fallback, str);
      FalsePositiveBenchmark(fallback, str);
    }

    ~LudoSeparator(){

    }


};

