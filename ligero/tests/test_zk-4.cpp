#include <cstdint>
#include <stdexcept>
#include <ctime>
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/fields/prime_base/fields_64.hpp>
#include "ligero/algebra/fft.hpp"
#include "ligero/algebra/polynomials/polynomial.hpp"
#include "ligero/protocols/ligero_iop.hpp"
#include <sys/time.h>

#include "ligero/tests/test_functions.hpp"

using namespace ligero;

int main()
{
    struct timeval start, end, total_start;
    typedef libff::Fields_64 FieldT;

    // 当type = 1时，是正确情况下的测试
    // 当type = 0时，是异常测试
    const std::size_t test_type = 1;

    // 修改input_size以完成输入规模扩展测试，可选参数为1, 2, 4, 8
    const std::size_t input_size = 8;

    const std::size_t RS_extra_dimensions = 2; // 2 or 3
    const std::size_t RS_col_extra_dimensions = 3;
    const std::vector<std::size_t> localization_parameter_array({1, 4, 4});
    const std::size_t security_parameter = 100; // P(soundness error)<2^(-security_parameter)
    const bool make_zk = true;
    const float height_width_ratio = 1 << 10;
    const std::size_t parallel_opt = 0; // 并行选项 0则不并行 其他为并行线程数

    gettimeofday(&total_start, nullptr);
    typedef sha256_two_to_one_hash_gadget<FieldT> HashT;
    protoboard<FieldT> pb;
    // #hash = tree_depth
    generate_merkle_check_read_example<FieldT, HashT>(pb, input_size);

    r1cs_constraint_system<FieldT> cs = pb.get_constraint_system();
    // EXPECT_TRUE(pb.is_satisfied());

    gettimeofday(&end, nullptr);
    const float circuit_generate = (end.tv_usec - total_start.tv_usec) / 1000000.0 + end.tv_sec - total_start.tv_sec;

    std::map<std::size_t, FieldT> public_input; // z
    public_input.insert(std::pair<std::size_t, FieldT>(pb.primary_input().size() + pb.auxiliary_input().size(), pb.auxiliary_input().back()));
    if (test_type == 0)
    {
        // 生成一个错误的输入
        public_input[0] = FieldT::random_element();
    }

    iop_protocol<FieldT> IOP;
    std::size_t num_variables = cs.num_variables();
    std::size_t num_constraints = cs.num_constraints();

    ligero_iop_parameters<FieldT> params(security_parameter,
                                         RS_extra_dimensions,
                                         RS_col_extra_dimensions,
                                         localization_parameter_array,
                                         height_width_ratio,
                                         make_zk,
                                         multiplicative_coset_type,
                                         num_constraints,
                                         num_variables,
                                         parallel_opt);

    ligero_iop<FieldT> proto(IOP,
                             cs,
                             params);

    gettimeofday(&end, nullptr);
    const float preprocessor = (end.tv_usec - total_start.tv_usec) / 1000000.0 + end.tv_sec - total_start.tv_sec;

    gettimeofday(&start, nullptr);
    proto.register_interactions();
    IOP.seal_interaction_registrations();
    proto.register_queries();
    IOP.seal_query_registrations();

    const float mt_arrange_time = proto.produce_oracle(pb.primary_input(), pb.auxiliary_input(), public_input);
    gettimeofday(&end, nullptr);
    const float oracle_generate = (end.tv_usec - start.tv_usec) / 1000000.0 + end.tv_sec - start.tv_sec - mt_arrange_time;

    gettimeofday(&start, nullptr);
    const float random_chuanzhi_time = proto.produce_proof(public_input);
    gettimeofday(&end, nullptr);
    const float proof_generate = (end.tv_usec - start.tv_usec) / 1000000.0 + end.tv_sec - start.tv_sec - random_chuanzhi_time;

    gettimeofday(&start, nullptr);
    bool correctness = proto.verifier_predicate(public_input);
    gettimeofday(&end, nullptr);
    //
    const float verify_time = (end.tv_usec - start.tv_usec) / 1000000.0 + end.tv_sec - start.tv_sec;
    const float total_time = (end.tv_usec - total_start.tv_usec) / 1000000.0 + end.tv_sec - total_start.tv_sec;
    if (correctness)
    {
        std::cout << "这是***基于安全多方计算的零知识证明模块-主要性能指标测试***" << std::endl;
        std::cout << "输入规模为深度为" << input_size << "的默克尔树！" << std::endl;
        // std::cout << "基于安全多方计算的零知识证明运行成功！" << std::endl;
        if (test_type == 1)
        {
            std::cout << "零知识证明验证通过！ " << std::endl;
            // std::cout
            //     << "基于安全多方计算的零知识证明总体功能测试通过！" << std::endl;
            // std::cout << "基于安全多方计算的零知识证明Ligerolight模块（模块1）执行通过！" << std::endl;
        }

        std::cout << "电路生成时间" << std::endl;
        std::cout << circuit_generate << " s" << std::endl;
        std::cout << "预处理时间" << std::endl;
        std::cout << preprocessor << " s" << std::endl;
        std::cout << "谕示生成时间" << std::endl;
        std::cout << oracle_generate << " s" << std::endl;
        std::cout << "证明生成时间时间" << std::endl;
        std::cout << proof_generate << " s" << std::endl;
        std::cout << "验证时间" << std::endl;
        std::cout << verify_time << "s" << std::endl;
        std::cout << "证明规模见控制台打印" << std::endl;
    }
    else
    {
        std::cout << "这是基于安全多方计算的零知识证明异常处理测试！" << std::endl;
        std::cout << "基于安全多方计算的零知识证明运行失败！" << std::endl;
        if (test_type == 0)
        {
            std::cout << "基于安全多方计算的零知识证明异常处理测试通过！" << std::endl;
        }
    }

    return 0;
}
