//#include <NormalEstimation.h>
//#include <ipsr_controller.h>
//
//template<typename REAL, int DIM>
//class IpsrNormalEstimator:public NormalEstimation<REAL,DIM>{
//private:
//    IPSR_Factory<REAL,DIM>* factory;
//    nlohmann::json config;
//
//public:
//    IpsrNormalEstimator(IPSR_Factory<REAL, DIM>* factory,int depth = 3,int iters = 10){
//        this->factory = new IPSR_Factory<REAL, DIM>(*factory);
//        this->factory->change_self_args("--depth",std::to_string(depth));
//        this->factory->change_self_args("--iters",std::to_string(iters));
//        config["depth"] = depth;
//        config["iters"] = iters;
//
//    }
//
//    int Estimate(ORIENTED_POINTS& op) {
//        IPSR_HANDLE_P hp = this->factory->create_ipsr_from_op(op, neumman_ind_ipsr, new RandomInit<REAL, DIM>(0));
//        int ep = -1;
//        do{
//            ep = hp->iter(); 
//        }while(ep != -1);
//        std::vector<int> tmp;
//        hp->update_target(op,hp->ixform.inverse(),tmp);
//        return 0;
//    }
//    nlohmann::json get_config() {
//        //config["ipsr_config"] = 
//        return config;
//    }
//};