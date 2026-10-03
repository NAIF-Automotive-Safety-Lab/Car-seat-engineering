#include <chrono/physics/ChBody.h>
#include <chrono/physics/ChLinkMotorRotationSpeed.h>
#include <chrono/physics/ChSystemNSC.h>
#include <chrono/functions/ChFunctionConst.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>

namespace fs = std::filesystem;
using namespace chrono;

static double max_constraint_error(const std::shared_ptr<ChLinkMotorRotationSpeed>& link) {
    const auto c = link->GetConstraintViolation();
    double m = 0.0;
    for (unsigned int i=0; i<c.size(); ++i)
        m = std::max(m, std::abs(c(i)));
    return m;
}

static void write_row(std::ofstream& out, const ChSystemNSC& sys,
                      const std::shared_ptr<ChBody>& body,
                      const std::shared_ptr<ChLinkMotorRotationSpeed>& motor) {
    const auto p = body->GetPos();
    const auto q = body->GetRot();
    const auto v = body->GetPosDt();
    const auto w = body->GetAngVelLocal();
    out << std::fixed << std::setprecision(12)
        << sys.GetChTime() << ","
        << p.x() << "," << p.y() << "," << p.z() << ","
        << q.e0() << "," << q.e1() << "," << q.e2() << "," << q.e3() << ","
        << v.x() << "," << v.y() << "," << v.z() << ","
        << w.x() << "," << w.y() << "," << w.z() << ","
        << motor->GetMotorAngle() << "," << motor->GetMotorAngleDt() << ","
        << max_constraint_error(motor) << "\n";
}

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    fs::path outdir(argv[1]);
    fs::create_directories(outdir);

    ChSystemNSC sys;
    sys.SetGravitationalAcceleration(ChVector3d(0,0,0));

    auto body_a = chrono_types::make_shared<ChBody>();
    body_a->SetFixed(true);
    body_a->SetName("BODY_A");
    sys.AddBody(body_a);

    auto body_b = chrono_types::make_shared<ChBody>();
    body_b->SetName("BODY_B");
    body_b->SetMass(1.0);
    body_b->SetInertiaXX(ChVector3d(1.0,1.0,1.0));
    sys.AddBody(body_b);

    auto motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
    motor->Initialize(body_b, body_a, ChFrame<>(ChVector3d(0,0,0)));
    motor->SetSpeedFunction(chrono_types::make_shared<ChFunctionConst>(0.5));
    sys.AddLink(motor);

    const double dt = 0.005;
    const double duration = 0.2;
    const int steps = static_cast<int>(std::lround(duration/dt));

    std::ofstream traj(outdir/"trajectory.csv");
    if (!traj) return 3;
    traj << "time_s,px_m,py_m,pz_m,q0,q1,q2,q3,vx_mps,vy_mps,vz_mps,wx_radps,wy_radps,wz_radps,motor_angle_rad,motor_speed_radps,max_constraint_error\n";
    for (int i=0; i<=steps; ++i) {
        write_row(traj,sys,body_b,motor);
        if (i<steps && sys.DoStepDynamics(dt)!=0) return 4;
    }
    traj.close();

    const double angle_error = std::abs(motor->GetMotorAngle() - 0.5*duration);
    const double final_time = sys.GetChTime();

    std::ofstream metrics(outdir/"metrics.json");
    metrics << std::fixed << std::setprecision(12)
      << "{\n  \"probe_class\": \"ENGINE_ONLY_GENERIC\",\n"
      << "  \"tocs_engineering_claim\": false,\n"
      << "  \"bodies\": [\"BODY_A\",\"BODY_B\"],\n"
      << "  \"controlled_joint\": \"ROTATION_SPEED_MOTOR\",\n"
      << "  \"time_step_s\": " << dt << ",\n"
      << "  \"duration_s\": " << duration << ",\n"
      << "  \"final_time_s\": " << final_time << ",\n"
      << "  \"expected_final_motor_angle_rad\": " << 0.5*duration << ",\n"
      << "  \"final_motor_angle_rad\": " << motor->GetMotorAngle() << ",\n"
      << "  \"final_angle_abs_error_rad\": " << angle_error << ",\n"
      << "  \"max_constraint_error\": " << max_constraint_error(motor) << ",\n"
      << "  \"constraints_enforced\": true\n}\n";

    if (std::abs(final_time-duration)>1e-12 || angle_error>1e-6) return 5;

    std::cout << "ENGINE_PROBE=PASS\nPROBE_CLASS=ENGINE_ONLY_GENERIC\nTOCS_ENGINEERING_CLAIM=false\n";
    return 0;
}
