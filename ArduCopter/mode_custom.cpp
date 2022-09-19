#include "Copter.h"
#include <AP_Motors/AP_MotorsMatrix.h>
#include <GCS_MAVLink/GCS.h>

// Function for hardcoding changes to MATLABs cntrl struct.
// Values can be accessed in the same fashion as in MATLAB, e.g.:
//     cntrl.sample_time = 42;
void ModeCustom::override_cntrl_params()
{

}

// init custom flight mode
bool ModeCustom::init(bool ignore_checks)
{
    override_cntrl_params();
    // initialize yaw to measured value
    const AP_AHRS_View &ahrs_ = attitude_control->get_ahrs();
    Quaternion attitude_vehicle_quat;
    ahrs_.get_quat_body_to_ned(attitude_vehicle_quat);
    yawInit = attitude_vehicle_quat.get_euler_yaw();
    updated_waypoints = true;
    // initialize position to measured value
    Vector3f position_NED;
    if (!ahrs_.get_relative_position_NED_home(position_NED)){
        position_NED[0] = 0;
        position_NED[1] = 0;
        position_NED[2] = 0; 
    }
    sInit[0] = position_NED[0];
    sInit[1] = position_NED[1];
    sInit[2] = position_NED[2];

    // tell the controller to use the initial conditions on the first time step
    custom_controller.initialize();

    return true;
}

// run custom flight mode
void ModeCustom::run()
{

    // Get stick inputs, -1 ... 1
    int16_t tr_max = 4500;
    // fetch roll and pitch inputs
    float roll_out_high = channel_roll->get_control_in();
    float roll_out = roll_out_high / tr_max;
    float pitch_out_high = channel_pitch->get_control_in();
    float pitch_out = pitch_out_high / tr_max;
    float throttle_control_high = channel_throttle->get_control_in();
    float throttle_control = throttle_control_high / 1000 * 2 - 1;
    // get pilot's desired yaw rate
    float yaw_out_high = channel_yaw->get_control_in();
    float yaw_out = yaw_out_high / tr_max;
    
    // Get measured values
    // Retrieve quaternion vehicle attitude
    const AP_AHRS_View& ahrs_ = attitude_control->get_ahrs();
    Quaternion attitude_vehicle_quat;
    ahrs_.get_quat_body_to_ned(attitude_vehicle_quat);
    // Get velocity relative to the ground in NED
    //bool check = ahrs_.have_inertial_nav(void);
    Vector3f velocity_NED;
    if (!ahrs_.get_velocity_NED(velocity_NED)) {
        velocity_NED[0] = 0;
        velocity_NED[1] = 0;
        velocity_NED[2] = 0;
    }

    Vector3f Omega_Kb_filt = ahrs_.get_gyro_latest();

    /* Info: gyro scaling is hard coded based on AP_InertialSensor::register_gyro in AP_InertialSensor.cpp.
    The scaling is applied in AP_InertialSensor_Backend::_notify_new_gyro_raw_sample in
    AP_InertialSensor_Backend.cpp. However, the variable gyro_filtered is overwritten during filtering.
    It seems that there is no non-filtered scaled angular velocity available as member variable.
    That is why the scaling is applied here.) */
    Vector3f Omega_Kb_raw = AP::ins().get_raw_gyro() / (INT16_MAX/radians(2000));

    float roll_angle = attitude_vehicle_quat.get_euler_roll();
    float pitch_angle = attitude_vehicle_quat.get_euler_pitch();
    float yaw_angle = attitude_vehicle_quat.get_euler_yaw();
    // Get position relative to the ground in NED
    Vector3f position_NED;
    if (!ahrs_.get_relative_position_NED_home(position_NED)) {
        position_NED[0] = 0;
        position_NED[1] = 0;
        position_NED[2] = 0;
    }


    // To do: spool states are currently based on copy from mode_stabilize
    if (!motors->armed()) {
        // Motors should be Stopped
        motors->set_desired_spool_state(AP_Motors::DesiredSpoolState::SHUT_DOWN);
    } else if (copter.ap.throttle_zero) {
        // Attempting to Land
        motors->set_desired_spool_state(AP_Motors::DesiredSpoolState::GROUND_IDLE);
    } else {
        motors->set_desired_spool_state(AP_Motors::DesiredSpoolState::THROTTLE_UNLIMITED);
    }
    // To do: spool states are currently based on copy from mode_stabilize
    switch (motors->get_spool_state()) {
    case AP_Motors::SpoolState::SHUT_DOWN:
        // Motors Stopped
        // To do!
        break;
    case AP_Motors::SpoolState::GROUND_IDLE:
        // Landed
        // To do!
        break;
    case AP_Motors::SpoolState::THROTTLE_UNLIMITED:
        // clear landing flag above zero throttle
        // To do!
    case AP_Motors::SpoolState::SPOOLING_UP:
    case AP_Motors::SpoolState::SPOOLING_DOWN:
        // do nothing
        // To do!
        break;
    }

    // assign commanded controller inputs to cmd struct
    ExtU rtU_;
    ExtY rtY_;

    rtU_.cmd.roll = roll_out;
    rtU_.cmd.pitch = pitch_out;
    rtU_.cmd.yaw = yaw_out;
    rtU_.cmd.thr = -throttle_control;
    rtU_.cmd.s_Kg_init[0] = sInit[0];
    rtU_.cmd.s_Kg_init[1] = sInit[1];
    rtU_.cmd.s_Kg_init[2] = sInit[2];
    rtU_.cmd.yaw_init = yawInit;
    for (int i=0;i<16;i++) {
        rtU_.cmd.RC_pwm[i] = g2.rc_channels.channel(i)->get_radio_in();
    }

    // assign or update waypoints
    // overwrite all custom controller waypoints with 5m above home position
    for (int k=0;k<max_num_of_matlab_waypoints;k++){
        rtU_.cmd.waypoints[4*k]   = 0.0f;
        rtU_.cmd.waypoints[4*k+1] = 0.0f;
        rtU_.cmd.waypoints[4*k+2] = -5.0f;
        rtU_.cmd.waypoints[4*k+3] = 0.0f;
    }
    int wp_count=0;
    // start with index j=1 because 1st Ardupilot waypoint is always home position
    for (int j=1;j<max_num_of_ardupilot_waypoints;j++){
        // assign only waypoints that are no "ghost waypoints", see declaration of waypoints
        if (abs(waypoints[j][0]) + abs(waypoints[j][1]) + abs(waypoints[j][2]) >= 0.01f){
            rtU_.cmd.waypoints[4*wp_count]   = waypoints[j][0]*0.01f; // convert cm to m
            rtU_.cmd.waypoints[4*wp_count+1] = waypoints[j][1]*0.01f; // convert cm to m
            rtU_.cmd.waypoints[4*wp_count+2] = waypoints[j][2]*0.01f; // convert cm to m
            rtU_.cmd.waypoints[4*wp_count+3] = waypoints[j][3]; // target velocity in m/s
            wp_count++;
        }
        if (wp_count>=max_num_of_matlab_waypoints){
            break;
        }
    }
    rtU_.cmd.num_waypoints = wp_count;
    rtU_.cmd.mission_change = updated_waypoints;
    updated_waypoints = false;


    // assign measured controller inputs to measure struct
    rtU_.measure.omega_Kb[0] = Omega_Kb_filt[0];
    rtU_.measure.omega_Kb[1] = Omega_Kb_filt[1];
    rtU_.measure.omega_Kb[2] = Omega_Kb_filt[2];
    rtU_.measure.Omega_Kb_raw[0] = Omega_Kb_raw[0];
    rtU_.measure.Omega_Kb_raw[1] = Omega_Kb_raw[1];
    rtU_.measure.Omega_Kb_raw[2] = Omega_Kb_raw[2];
    rtU_.measure.q_bg[0] = attitude_vehicle_quat.q1;
    rtU_.measure.q_bg[1] = attitude_vehicle_quat.q2;
    rtU_.measure.q_bg[2] = attitude_vehicle_quat.q3;
    rtU_.measure.q_bg[3] = attitude_vehicle_quat.q4;
    rtU_.measure.EulerAngles[0] = roll_angle;
    rtU_.measure.EulerAngles[1] = pitch_angle;
    rtU_.measure.EulerAngles[2] = yaw_angle;
    rtU_.measure.a_Kg[0] = ahrs_.get_accel_ef_blended().x;
    rtU_.measure.a_Kg[1] = ahrs_.get_accel_ef_blended().y;
    rtU_.measure.a_Kg[2] = ahrs_.get_accel_ef_blended().z;
    rtU_.measure.V_Kg[0] = velocity_NED[0];
    rtU_.measure.V_Kg[1] = velocity_NED[1];
    rtU_.measure.V_Kg[2] = velocity_NED[2];
    rtU_.measure.s_Kg[0] = position_NED[0];
    rtU_.measure.s_Kg[1] = position_NED[1];
    rtU_.measure.s_Kg[2] = position_NED[2];
    rtU_.measure.lla[0] = copter.current_loc.lat;
    rtU_.measure.lla[1] = copter.current_loc.lng;
    rtU_.measure.lla[2] = copter.current_loc.alt;
    
    // run Simulink controller
    custom_controller.rtU = rtU_;
    custom_controller.step();
    rtY_ = custom_controller.rtY;

    // DEBUGGING:
    // Send measure bus to Simulink (uncomment line 3 in mode.h)
    #ifdef Custom_Matlab_Output
        socket_debug.sendto(&rtU_.measure, sizeof(rtU_.measure), _debug_address, _debug_port); 
    #endif

    // log data
    AP::logger().Write(
        "ML", "TimeUS,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15",
        "Qfffffffffffffff",
        AP_HAL::micros64(),
        (double)custom_controller.rtY.logs[0],
        (double)custom_controller.rtY.logs[1],
        (double)custom_controller.rtY.logs[2],
        (double)custom_controller.rtY.logs[3],
        (double)custom_controller.rtY.logs[4],
        (double)custom_controller.rtY.logs[5],
        (double)custom_controller.rtY.logs[6],
        (double)custom_controller.rtY.logs[7],
        (double)custom_controller.rtY.logs[8],
        (double)custom_controller.rtY.logs[9],
        (double)custom_controller.rtY.logs[10],
        (double)custom_controller.rtY.logs[11],
        (double)custom_controller.rtY.logs[12],
        (double)custom_controller.rtY.logs[13],
        (double)custom_controller.rtY.logs[14]);

    // set outputs in the same order as Simulink
    for (int i=0;i<8;i++) {
        motors->set_custom_input( i, rtY_.u[i] );
    }
}


void ModeCustom::add_waypoint(uint16_T index,Vector3f location){
        waypoints[index][0] = location.x;
        waypoints[index][1] = location.y;
        waypoints[index][2] = -location.z;
        waypoints[index][3] = 0.0f;
}

void ModeCustom::add_speed(uint16_T index, float V_k){
    if(abs(waypoints[index-1][0] + waypoints[index-1][1] + waypoints[index-1][2]) >= 0.1f){
        waypoints[index][0] = 0.0f;
        waypoints[index][1] = 0.0f;
        waypoints[index][2]     = 0.0f;
        waypoints[index-1][3] = V_k;
    }

}