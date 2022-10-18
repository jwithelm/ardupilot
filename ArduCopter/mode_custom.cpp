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

    // init custom logging
    log_setup(log_config);

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
    Vector3f position_NED_origin;
    if (!ahrs_.get_relative_position_NED_origin(position_NED_origin)) {
        position_NED_origin[0] = 0;
        position_NED_origin[1] = 0;
        position_NED_origin[2] = 0;
    }
    float voltage = copter.battery.voltage();


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
    rtU_.measure.s_Kg_origin[0] = position_NED[0];
    rtU_.measure.s_Kg_origin[1] = position_NED_origin[1];
    rtU_.measure.s_Kg_origin[2] = position_NED_origin[2];
    rtU_.measure.lla[0] = copter.current_loc.lat;
    rtU_.measure.lla[1] = copter.current_loc.lng;
    rtU_.measure.lla[2] = copter.current_loc.alt;
    rtU_.measure.V_bat = voltage;
    
    // run Simulink controller
    custom_controller.rtU = rtU_;
    custom_controller.step();
    rtY_ = custom_controller.rtY;

    // DEBUGGING:
    // Send measure bus to Simulink (uncomment line 3 in mode.h)
    #ifdef Custom_Matlab_Output
        socket_debug.sendto(&rtU_.measure, sizeof(rtU_.measure), _debug_address, _debug_port); 
    #endif

    // log signals
    for (int i=0;i<num_log_batches;i++) {
        char label[label_length[i]+1];
        get_log_label(i+1, label);
        char batch_name[batch_name_length[i]+1];
        get_log_batch_name(i+1, batch_name);
        write_log_custom(batch_name, label,
            &custom_controller.rtY.logs[log_signal_idx_cumsum[i]],
            log_config[i].num_signals);
    }

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
        waypoints[index][2] = 0.0f;
        waypoints[index-1][3] = V_k;
    }
}

void ModeCustom::log_setup(const logConfigBus log_config_in[]) {
    set_log_batch_names(log_config_in);
    set_log_labels(log_config_in);
    set_log_signal_idx_cumsum(log_config_in);
};

void ModeCustom::set_log_batch_names(const logConfigBus log_config_in[]) {
    for (int i=0;i<num_log_batches;i++) {
        memcpy(&(batch_name_full[i][0]), &(log_config_in[i].batch_name), max_batch_name_length);
        batch_name_length[i]=max_batch_name_length;
        for (int j=0;j<max_batch_name_length;j++) {
            if (batch_name_full[i][j] == 1) {
                batch_name_length[i] --;
            }
        }
    }
};

void ModeCustom::set_log_labels(const logConfigBus log_config_in[]){
    signal_name_t current_name_int;
    int signal_name_length;
    for (int i=0;i<num_log_batches;i++) {
        const uint8_t *log_names_in = log_config_in[i].signal_names;
        memcpy(&(label_full[i][0]), &("TimeUS"), 6);
        label_length[i]=6;
        for (int j=0;j<log_config_in[i].num_signals;j++) {
            signal_name_length = max_signal_name_length;
            memcpy(&(label_full[i][label_length[i]]),&(","),1);
            label_length[i] ++;
            extract_one_signal_name(log_names_in, j+1, current_name_int);
            for (int k=0;k<max_signal_name_length;k++) {
                if (current_name_int[k]==1) {
                    signal_name_length --;
                }
            }
            memcpy(&(label_full[i][label_length[i]]),&current_name_int,signal_name_length);
            label_length[i] += signal_name_length;
        }
    }
};

void ModeCustom::set_log_signal_idx_cumsum(const logConfigBus log_config_in[]){
    log_signal_idx_cumsum[0] = 0;
    for (int i=1;i<num_log_batches;i++) {
        log_signal_idx_cumsum[i] = log_signal_idx_cumsum[i-1] + log_config_in[i-1].num_signals;
    }
};

void ModeCustom::write_log_custom(const char *name, const char *labels, float *sf, int size) {
    double s[size];
    for (int i=0; i<size; i++) {
        s[i] = (double)sf[i];
    }
    if (size==0) {
        return;
    } else if (size==1) {
        AP::logger().Write(name, labels,"Qf",AP_HAL::micros64(),s[0]);
    } else if (size==2) {
        AP::logger().Write(name, labels,"Qff",AP_HAL::micros64(),s[0],s[1]);
    } else if (size==3) {
        AP::logger().Write(name, labels,"Qfff",AP_HAL::micros64(),s[0],s[1],s[2]);
    } else if (size==4) {
        AP::logger().Write(name, labels,"Qffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3]);
    } else if (size==5) {
        AP::logger().Write(name, labels,"Qfffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4]);
    } else if (size==6) {
        AP::logger().Write(name, labels,"Qffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5]);
    } else if (size==7) {
        AP::logger().Write(name, labels,"Qfffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6]);
    } else if (size==8) {
        AP::logger().Write(name, labels,"Qffffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7]);
    } else if (size==9) {
        AP::logger().Write(name, labels,"Qfffffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8]);
    } else if (size==10) {
        AP::logger().Write(name, labels,"Qffffffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9]);
    } else if (size==11) {
        AP::logger().Write(name, labels,"Qfffffffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10]);
    } else if (size==12) {
        AP::logger().Write(name, labels,"Qffffffffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11]);
    } else if (size==13) {
        AP::logger().Write(name, labels,"Qfffffffffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12]);
    } else if (size==14) {
        AP::logger().Write(name, labels,"Qffffffffffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13]);
    } else if (size==15) {
        AP::logger().Write(name, labels,"Qfffffffffffffff",AP_HAL::micros64(),s[0],s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9],s[10],s[11],s[12],s[13],s[14]);
    }
};

void ModeCustom::extract_one_signal_name(const uint8_t log_names_int[], int number, signal_name_t &log_name){
    int idx = (number-1)*max_signal_name_length;
    for (int i=0;i<max_signal_name_length;i++) {
        log_name[i] = log_names_int[idx+i];
    }
};

void ModeCustom::get_log_label(int batch_number, char *label) {
    memcpy(label,&(label_full[batch_number-1]),label_length[batch_number-1]+1);
    label[label_length[batch_number-1]]=0;
};

void ModeCustom::get_log_batch_name(int batch_number, char *name) {
    memcpy(name,&(batch_name_full[batch_number-1]),batch_name_length[batch_number-1]+1);
    name[batch_name_length[batch_number-1]]=0;
};
