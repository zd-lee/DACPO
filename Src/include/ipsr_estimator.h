//
//private:
//    IPSR_Factory<REAL,DIM>* factory;
//    nlohmann::json config;
//
//public:
//        this->factory = new IPSR_Factory<REAL, DIM>(*factory);
//        this->factory->change_self_args("--depth",std::to_string(depth));
//        this->factory->change_self_args("--iters",std::to_string(iters));
//        config["depth"] = depth;
//        config["iters"] = iters;
//
//
//        IPSR_HANDLE_P hp = this->factory->create_ipsr_from_op(op, neumman_ind_ipsr, new RandomInit<REAL, DIM>(0));
//        }while(ep != -1);
//        hp->update_target(op,hp->ixform.inverse(),tmp);
//    nlohmann::json get_config() {
//        //config["ipsr_config"] = 
//};