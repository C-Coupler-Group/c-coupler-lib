/***************************************************************
  *  Copyright (c) 2017, Tsinghua University.
  *  This is a source file of C-Coupler.
  *  This file was initially finished by Dr. Li Liu. 
  *  If you have any problem, 
  *  please contact Dr. Li Liu via liuli-cess@tsinghua.edu.cn
  ***************************************************************/


#ifndef COUPLING_TIME_MGT
#define COUPLING_TIME_MGT


#define FREQUENCY_UNIT_STEP             "step"
#define FREQUENCY_UNIT_SECOND           "second"
#define FREQUENCY_UNIT_MINUTE           "minute"
#define FREQUENCY_UNIT_HOUR             "hour"
#define FREQUENCY_UNIT_DAY              "day"
#define FREQUENCY_UNIT_MONTH            "month"
#define FREQUENCY_UNIT_YEAR             "year"
#define FREQUENCY_UNIT_STEPS            "steps"
#define FREQUENCY_UNIT_SECONDS          "seconds"
#define FREQUENCY_UNIT_MINUTES          "minutes"
#define FREQUENCY_UNIT_HOURS            "hours"
#define FREQUENCY_UNIT_DAYS             "days"
#define FREQUENCY_UNIT_MONTHS           "months"
#define FREQUENCY_UNIT_YEARS            "years"
#define FREQUENCY_UNIT_NSTEP            "nstep"
#define FREQUENCY_UNIT_NSECOND          "nsecond"
#define FREQUENCY_UNIT_NMINUTE          "nminute"
#define FREQUENCY_UNIT_NHOUR            "nhour"
#define FREQUENCY_UNIT_NDAY             "nday"
#define FREQUENCY_UNIT_NMONTH           "nmonth"
#define FREQUENCY_UNIT_NYEAR            "nyear"
#define FREQUENCY_UNIT_NSTEPS           "nsteps"
#define FREQUENCY_UNIT_NSECONDS         "nseconds"
#define FREQUENCY_UNIT_NMINUTES         "nminutes"
#define FREQUENCY_UNIT_NHOURS           "nhours"
#define FREQUENCY_UNIT_NDAYS            "ndays"
#define FREQUENCY_UNIT_NMONTHS          "nmonths"
#define FREQUENCY_UNIT_NYEARS           "nyears"

#define RUNTYPE_INITIAL                 "initial"
#define RUNTYPE_CONTINUE                "continue"
#define RUNTYPE_BRANCH                  "branch"
#define RUNTYPE_HYBRID                  "hybrid"

#define SECONDS_PER_DAY                 ((long)86400)
#define NUM_MONTH_PER_YEAR              ((long)12)
#define NUM_DAYS_PER_NONLEAP_YEAR       ((long)365)
#define NUM_DAYS_PER_LEAP_YEAR          ((long)366)


#define IS_TIME_UNIT_STEP(unit)         (words_are_the_same(unit,FREQUENCY_UNIT_STEP) || words_are_the_same(unit,FREQUENCY_UNIT_STEPS) || words_are_the_same(unit,FREQUENCY_UNIT_NSTEP) || words_are_the_same(unit,FREQUENCY_UNIT_NSTEPS))
#define IS_TIME_UNIT_SECOND(unit)       (words_are_the_same(unit,FREQUENCY_UNIT_SECOND) || words_are_the_same(unit,FREQUENCY_UNIT_SECONDS) || words_are_the_same(unit,FREQUENCY_UNIT_NSECOND) || words_are_the_same(unit,FREQUENCY_UNIT_NSECONDS))
#define IS_TIME_UNIT_MINUTE(unit)       (words_are_the_same(unit,FREQUENCY_UNIT_MINUTE) || words_are_the_same(unit,FREQUENCY_UNIT_MINUTES) || words_are_the_same(unit,FREQUENCY_UNIT_NMINUTE) || words_are_the_same(unit,FREQUENCY_UNIT_NMINUTES))
#define IS_TIME_UNIT_HOUR(unit)         (words_are_the_same(unit,FREQUENCY_UNIT_HOUR) || words_are_the_same(unit,FREQUENCY_UNIT_HOURS) || words_are_the_same(unit,FREQUENCY_UNIT_NHOUR) || words_are_the_same(unit,FREQUENCY_UNIT_NHOURS))
#define IS_TIME_UNIT_DAY(unit)          (words_are_the_same(unit,FREQUENCY_UNIT_DAY) || words_are_the_same(unit,FREQUENCY_UNIT_DAYS) || words_are_the_same(unit,FREQUENCY_UNIT_NDAY) || words_are_the_same(unit,FREQUENCY_UNIT_NDAYS))
#define IS_TIME_UNIT_MONTH(unit)        (words_are_the_same(unit,FREQUENCY_UNIT_MONTH) || words_are_the_same(unit,FREQUENCY_UNIT_MONTHS) || words_are_the_same(unit,FREQUENCY_UNIT_NMONTH) || words_are_the_same(unit,FREQUENCY_UNIT_NMONTHS))
#define IS_TIME_UNIT_YEAR(unit)         (words_are_the_same(unit,FREQUENCY_UNIT_YEAR) || words_are_the_same(unit,FREQUENCY_UNIT_YEARS) || words_are_the_same(unit,FREQUENCY_UNIT_NYEAR) || words_are_the_same(unit,FREQUENCY_UNIT_NYEARS))


#define TIME_UNIT_STRING_STEP           "\"step\", \"steps\", \"nstep\", \"nsteps\""
#define TIME_UNIT_STRING_SECOND         "\"second\", \"seconds\", \"nsecond\", \"nseconds\""
#define TIME_UNIT_STRING_MINUTE         "\"minute\", \"minutes\", \"nminute\", \"nminutes\""
#define TIME_UNIT_STRING_HOUR           "\"hour\", \"hours\", \"nhour\", \"nhours\""
#define TIME_UNIT_STRING_DAY            "\"day\", \"days\", \"nday\", \"ndays\""
#define TIME_UNIT_STRING_MONTH          "\"month\", \"months\", \"nmonth\", \"nmonths\""
#define TIME_UNIT_STRING_YEAR           "\"year\", \"years\", \"nyear\", \"nyears\""

#define RUNTYPE_MARK_INITIAL            1
#define RUNTYPE_MARK_CONTINUE           2
#define RUNTYPE_MARK_BRANCH             3
#define RUNTYPE_MARK_HYBRID             4


#include "common_utils.h"
#include <vector>


class Time_mgt;


class Coupling_timer
{
    private:
        friend class Time_mgt;
        char frequency_unit[NAME_STR_SIZE];
        int frequency_count;
        int local_lag_count;
        int remote_lag_count;
        int timer_id;
        int comp_id;
        std::vector<Coupling_timer*> children;
        Time_mgt *comp_time_mgr;
        int or_or_and;
        
    public:
        Coupling_timer(int, int, int*, int, int, const char *);
        Coupling_timer(int, int, const char*, int, int, int, const char*);
        Coupling_timer(int, int, Coupling_timer*);
        Coupling_timer(const char*, long &, int, bool, bool &);
        ~Coupling_timer();
        bool is_timer_on();
        bool is_timer_on(int, int, int, int, int, int, int, int, int, int);
        int get_timer_id() { return timer_id; }
        int get_comp_id() { return comp_id; }
        int get_frequency_count() { return frequency_count; }
        int get_local_lag_count() { return local_lag_count; } 
        int get_remote_lag_count() { return remote_lag_count; }
        const char *get_frequency_unit() { return frequency_unit; }
        void write_timer_into_array(char **, long &, long &);
        void get_time_of_next_timer_on(Time_mgt *, int, int, int, int, int, int, int &, int &, int &, bool);
        void reset_remote_lag_count() { remote_lag_count = 0; }
        void check_timer_format();
        bool is_the_same_with(Coupling_timer *);
};


class Timer_mgt
{
    private:
        std::vector<Coupling_timer*> timers;

    public:
        Timer_mgt() {};
        ~Timer_mgt();
        bool check_is_legal_timer_id(int);
        Coupling_timer *get_timer(int);
        int define_timer(int, const char*, int, int, int, const char*);
        int define_timer(int, int*, int, int, int, const char*);
        int define_timer(int, Coupling_timer*);
        bool is_timer_on(int, const char*);
        void add_timer(Coupling_timer *);
};


class Time_mgt
{
    private:
        int start_year;
        int start_month;
        int start_day;
        int start_second;
        int restart_second;
        int previous_year;   
        int previous_month;
        int previous_day; 
        int previous_second; 
        int current_year;   
        int current_month;
        int current_day; 
        int current_second;
        int reference_year;
        int reference_month;
        int reference_day;
        int stop_year;
        int stop_month;
        int stop_day;
        int stop_second;
        int time_step_in_second; 
        int current_num_elapsed_day;
        int start_num_elapsed_day;
        int restart_num_elapsed_day;        
        long common_restart_full_time;
        long restart_full_time;
        int stop_num_elapsed_day;
        int current_step_id;
        int restarted_step_id;
        long num_total_steps;
        bool leap_year_on;
        Coupling_timer *restart_timer;
        int comp_id;
        char case_name[NAME_STR_SIZE];
        char exp_model_name[NAME_STR_SIZE];
        char case_desc[NAME_STR_SIZE];
        char run_type[NAME_STR_SIZE];
        int runtype_mark;
        char stop_option[NAME_STR_SIZE];
        char rest_freq_unit[NAME_STR_SIZE];
        int rest_freq_count;
        char rest_refcase[NAME_STR_SIZE];
        int rest_refdate;
        int rest_refsecond;
        bool advance_time_synchronized;
        bool time_has_been_advanced;
        int stop_n;

    public:
        Time_mgt() {}
        Time_mgt(int, const char *, bool);
        ~Time_mgt();
        void initialize_to_start_time();
        void advance_model_time(const char*, bool);
        void advance_time(int &, int &, int &, int &, int &, int);
        int get_current_year() { return current_year; }
        int get_current_month() { return current_month; }
        int get_current_day() { return current_day; }
        int get_current_hour() { return current_second /3600; }
        int get_current_minute() { return (current_second % 3600) / 60; }
        int get_current_second() { return current_second; }
        int get_time_step_in_second() { return time_step_in_second; }
        int get_stop_year() { return stop_year; }
        int get_stop_month() { return stop_month; }
        int get_stop_day() { return stop_day; }
        int get_stop_second() { return stop_second; }
        int get_start_year() { return start_year; }
        int get_start_month() { return start_month; }
        int get_start_day() { return start_day; }
        int get_start_second() { return start_second; }
        int get_start_num_elapsed_day() { return start_num_elapsed_day; }
        int get_stop_num_elapsed_day() { return stop_num_elapsed_day; }
        bool is_timer_on(const char *, int, int);
        bool check_is_model_run_finished();
        bool check_is_coupled_run_restart_time();
        double get_double_current_calendar_time(int, const char*);
        float get_float_current_calendar_time(int, const char*);
        long get_start_full_time();
        long get_previous_full_time();
        long get_current_full_time();
        int get_current_date();
        int get_current_num_time_step();
        long get_num_total_step() { return num_total_steps; }
        int get_current_num_days_in_year();
        void check_timer_format(const char*, int, int, int, bool, const char*);
        bool check_time_consistency_between_components(long);
        long calculate_elapsed_day(int, int, int);
		long get_elapsed_day_from_full_time(long);
        void get_elapsed_days_from_start_date(int*, int*);
        void get_elapsed_days_from_reference_date(int*, int*);
        void get_current_time(int&, int&, int&, int&, int, const char*);
        bool check_is_time_legal(int, int, int, int, const char*);
        bool get_is_leap_year_on() { return leap_year_on; }
        int get_comp_id() { return comp_id; }
        Time_mgt *clone_time_mgr(int);
        bool set_time_step_in_second(int, const char*, bool);
        bool is_a_leap_year(int);
        void build_restart_timer();
        int get_current_step_id() { return current_step_id; }
        void check_consistency_of_current_time(int, int, const char*);
        int get_current_num_elapsed_day() { return current_num_elapsed_day; }
        bool is_time_out_of_execution(long);
        void write_time_mgt_into_array(char **, long &, long &);
        void import_restart_data(const char*, long&, const char *, bool);
        bool is_restart_timer_on();
        const char *get_case_name() { return case_name; }
        const char *get_run_type() { return run_type; }
        const char *get_rest_refcase() { return rest_refcase; }
        int get_rest_refdate() { return rest_refdate; }
        int get_rest_refsecond() { return rest_refsecond; }
        bool is_first_restart_step() { return current_step_id == restarted_step_id; }
        void calculate_stop_time(int, int, int, int);
        bool get_time_has_been_advanced() { return time_has_been_advanced; }
        void reset_current_time_to_start_time(const char*);
        long get_restart_full_time() { return restart_full_time; }
        long get_common_restart_full_time() { return common_restart_full_time; }
        int get_runtype_mark() { return runtype_mark; }
        long determine_continue_run_restart_time();
};


class Components_time_mgt
{
    private:
        std::vector<Time_mgt*> components_time_mgrs;

    public:
        Components_time_mgt() {}
        ~Components_time_mgt();
        Time_mgt *get_time_mgr(int);
        void define_root_comp_time_mgr(int, const char*);
        void clone_parent_comp_time_mgr(int, int, const char*);
        void set_component_time_step(int, int, const char*);
        void advance_component_time(int, const char*);
        bool is_model_run_ended(int, const char*);
        void check_component_current_time(int, int, int, const char*);
};


extern int num_days_of_month_of_nonleap_year[];
extern int num_days_of_month_of_leap_year[];


#endif
