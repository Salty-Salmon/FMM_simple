#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <chrono>

class Stopwatch{
private:
    struct Note{
        std::string name;
        unsigned int counter;
        unsigned long long duration;

        Note(std::string& name, unsigned long long duration);
    };
    std::vector<Note> notes;
public:
    static unsigned long long get_nanoseconds();
    void propose_note(std::string& name, unsigned long long duration);
    friend std::ostream& operator<<(std::ostream& out, Stopwatch& stopwatch);
};
extern Stopwatch stopwatch;

class Stopwatch_call{
private:
    std::string name;
    unsigned long long t_0;
public:
    static Stopwatch *stopwatch_p;
    Stopwatch_call(std::string name);
    ~Stopwatch_call();
};
