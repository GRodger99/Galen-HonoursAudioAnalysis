#include <vamp-hostsdk/PluginHostAdapter.h>
#include <vamp-hostsdk/PluginInputDomainAdapter.h>
#include <vamp-hostsdk/PluginLoader.h>
#include <vamp-hostsdk/vamp-hostsdk.h>
#include <iostream>

#include <fstream>
#include <set>
#include <sndfile.h>

#include <cstring>
#include <cstdlib>
#include <string>
#include <cmath>
#include <ctime>
#include <chrono>
#include <algorithm>

#include <filesystem>
#include <Windows.h>

using namespace std;
namespace fs = std::filesystem;
using Vamp::Plugin;
using Vamp::PluginHostAdapter;
using Vamp::RealTime;
using Vamp::HostExt::PluginLoader;
using Vamp::HostExt::PluginWrapper;
using Vamp::HostExt::PluginInputDomainAdapter;

float lastDuration;
float songDuration;
struct SegmentDetails
{
	
	bool repeats = false;
	bool HighIntensity = false;
	float AverageIntensity = 0;
};
typedef std::pair<int, SegmentDetails> int_Seg;


#define HOST_VERSION "1"
class Timer				//Timer class taken from this github repo.  https://gist.github.com/gongzhitaao/7062087
{
public:
	Timer() : beg_(clock_::now()) {}
	void reset() { beg_ = clock_::now(); }
	double elapsed() const {
		return std::chrono::duration_cast<second_>
			(clock_::now() - beg_).count();
	}

private:
	typedef std::chrono::high_resolution_clock clock_;
	typedef std::chrono::duration<double, std::ratio<1> > second_;
	std::chrono::time_point<clock_> beg_;
};

std::vector<float> combined;
std::vector<float> miliSeconds;
std::vector<int> seconds;

std::vector<float> combinedSegment;
std::vector<float> miliSecondsSegment;
std::vector<int> secondsSegment;
std::vector<float> segmentNames;

std::vector<float> combinedIntensity;
std::vector<float> miliSecondsIntensity;
std::vector<int> secondsIntensity;
std::vector<float> IntensityValues;

bool Beats = true;
bool Segments = false;
bool intensity = false;
enum Verbosity {
	PluginIds,
	PluginOutputIds,
	PluginInformation,
	PluginInformationDetailed
};
#define max(a,b)            (((a) > (b)) ? (a) : (b))

static double
toSeconds(const RealTime& time)
{
	return time.sec + double(time.nsec + 1) / 1000000000.0;
}
void
printFeatures(int frame, int sr,
	const Plugin::OutputDescriptor& output, int outputNo,
	const Plugin::FeatureSet& features, ofstream* out, bool useFrames)
{
	static int featureCount = -1;

	if (features.find(outputNo) == features.end()) return;

	for (size_t i = 0; i < features.at(outputNo).size(); ++i) {

		const Plugin::Feature& f = features.at(outputNo).at(i);

		bool haveRt = false;
		RealTime rt;

		if (output.sampleType == Plugin::OutputDescriptor::VariableSampleRate) {
			rt = f.timestamp;
			haveRt = true;
		}
		else if (output.sampleType == Plugin::OutputDescriptor::FixedSampleRate) {
			int n = featureCount + 1;
			if (f.hasTimestamp) {
				n = int(round(toSeconds(f.timestamp) * output.sampleRate));
			}
			rt = RealTime::fromSeconds(double(n) / output.sampleRate);
			haveRt = true;
			featureCount = n;
		}

		if (useFrames) {

			int displayFrame = frame;

			if (haveRt) {
				displayFrame = RealTime::realTime2Frame(rt, sr);
			}

			(out ? *out : cout) << displayFrame;

			if (f.hasDuration) {
				displayFrame = RealTime::realTime2Frame(f.duration, sr);
				(out ? *out : cout) << "," << displayFrame;
			}

			(out ? *out : cout) << ":";

		}
		else {

			if (!haveRt) {
				rt = RealTime::frame2RealTime(frame, sr);
			}

			(out ? *out : cout) << rt.toString();
			bool pad = false;

			/////////////////////////////////////////////////////////////////////////
			if (out && Beats == true)
			{
				seconds.push_back(rt.sec);
				miliSeconds.push_back(rt.msec());

			}
			if (out && Segments == true)
			{
				secondsSegment.push_back(rt.sec);
				miliSecondsSegment.push_back(rt.msec());

			}
			if (out && intensity == true)
			{
				secondsIntensity.push_back(rt.sec);
				miliSecondsIntensity.push_back(rt.msec());
			}

			//////////////////////////////////////////////////////////////////////////////
			if (f.hasDuration) {
				rt = f.duration;
				(out ? *out : cout) << "," << rt.toString();
				lastDuration = rt.sec;
			}

			(out ? *out : cout) << ":";
		}
	if (out && Segments == true)
	{
		for (unsigned int j = 0; j < f.values.size(); ++j) {
			(out ? *out : cout) << " " << f.values[j];
			segmentNames.push_back(f.values[j]);
		}
		(out ? *out : cout) << " " << f.label;

		(out ? *out : cout) << endl;
	}
	else if(out && intensity == true)
	{
		for (unsigned int j = 0; j < f.values.size(); ++j) {
			(out ? *out : cout) << " " << f.values[j];
			IntensityValues.push_back(f.values[j]);
		}
		(out ? *out : cout) << " " << f.label;

		(out ? *out : cout) << endl;
	}
	else 
	{
		for (unsigned int j = 0; j < f.values.size(); ++j) {
			(out ? *out : cout) << " " << f.values[j];
		}
		(out ? *out : cout) << " " << f.label;

		(out ? *out : cout) << endl;
	}
}
}




int run(string myname, string soname, string id,
	string output, int outputNo, string wavname,
	string outfilename, bool useFrames)
{


	{
		PluginLoader* loader = PluginLoader::getInstance();
	
		PluginLoader::PluginKey key = loader->composePluginKey(soname, id);

		SNDFILE* sndfile;
		SF_INFO sfinfo;
		
		memset(&sfinfo, 0, sizeof(SF_INFO));

		sndfile = sf_open(wavname.c_str(), SFM_READ, &sfinfo);
		if (!sndfile) {
			cerr << myname << ": ERROR: Failed to open input file \""
				<< wavname << "\": " << sf_strerror(sndfile) << endl;
			return 1;
		}

		ofstream* out = 0;
		if (outfilename != "") {
			out = new ofstream(outfilename.c_str(), ios::out);
			if (!*out) {
				cerr << myname << ": ERROR: Failed to open output file \""
					<< outfilename << "\" for writing" << endl;
				delete out;
				return 1;
			}
		}
		
		Plugin* plugin = loader->loadPlugin
		(key, sfinfo.samplerate, PluginLoader::ADAPT_ALL_SAFE);
		if (!plugin) {
			cerr << myname << ": ERROR: Failed to load plugin \"" << id
				<< "\" from library \"" << soname << "\"" << endl;
			sf_close(sndfile);
			if (out) {
				out->close();
				delete out;
			}
			return 1;
		}

		cerr << "Running plugin: \"" << plugin->getIdentifier() << "\"..." << endl;

		// Note that the following would be much simpler if we used a
		// PluginBufferingAdapter as well -- i.e. if we had passed
		// PluginLoader::ADAPT_ALL to loader->loadPlugin() above, instead
		// of ADAPT_ALL_SAFE.  Then we could simply specify our own block
		// size, keep the step size equal to the block size, and ignore
		// the plugin's bleatings.  However, there are some issues with
		// using a PluginBufferingAdapter that make the results sometimes
		// technically different from (if effectively the same as) the
		// un-adapted plugin, so we aren't doing that here.  See the
		// PluginBufferingAdapter documentation for details.

		int blockSize = plugin->getPreferredBlockSize();
		int stepSize = plugin->getPreferredStepSize();

		if (blockSize == 0) {
			blockSize = 1024;
		}
		if (stepSize == 0) {
			if (plugin->getInputDomain() == Plugin::FrequencyDomain) {
				stepSize = blockSize / 2;
			}
			else {
				stepSize = blockSize;
			}
		}
		else if (stepSize > blockSize) {
			cerr << "WARNING: stepSize " << stepSize << " > blockSize " << blockSize << ", resetting blockSize to ";
			if (plugin->getInputDomain() == Plugin::FrequencyDomain) {
				blockSize = stepSize * 2;
			}
			else {
				blockSize = stepSize;
			}
			cerr << blockSize << endl;
		}
		int overlapSize = blockSize - stepSize;
		sf_count_t currentStep = 0;
		int finalStepsRemaining = max(1, (blockSize / stepSize) - 1); // at end of file, this many part-silent frames needed after we hit EOF

		int channels = sfinfo.channels;

		float* filebuf = new float[blockSize * channels];
		float** plugbuf = new float* [channels];
		for (int c = 0; c < channels; ++c) plugbuf[c] = new float[blockSize + 2];

		cerr << "Using block size = " << blockSize << ", step size = "
			<< stepSize << endl;

		// The channel queries here are for informational purposes only --
		// a PluginChannelAdapter is being used automatically behind the
		// scenes, and it will take case of any channel mismatch

		int minch = plugin->getMinChannelCount();
		int maxch = plugin->getMaxChannelCount();
		cerr << "Plugin accepts " << minch << " -> " << maxch << " channel(s)" << endl;
		cerr << "Sound file has " << channels << " (will mix/augment if necessary)" << endl;

		Plugin::OutputList outputs = plugin->getOutputDescriptors();
		Plugin::OutputDescriptor od;
		Plugin::FeatureSet features;
		Plugin::ParameterList parameters = plugin->getParameterDescriptors();
		int returnValue = 1;
		int progress = 0;

		RealTime rt;
		PluginWrapper* wrapper = 0;
		RealTime adjustment = RealTime::zeroTime;

		if (outputs.empty()) {
			cerr << "ERROR: Plugin has no outputs!" << endl;
			goto done;
		}

		if (outputNo < 0) {

			for (size_t oi = 0; oi < outputs.size(); ++oi) {
				if (outputs[oi].identifier == output) {
					outputNo = oi;
					break;
				}
			}

			if (outputNo < 0) {
				cerr << "ERROR: Non-existent output \"" << output << "\" requested" << endl;
				goto done;
			}

		}
		else {

			if (int(outputs.size()) <= outputNo) {
				cerr << "ERROR: Output " << outputNo << " requested, but plugin has only " << outputs.size() << " output(s)" << endl;
				goto done;
			}
		}

		od = outputs[outputNo];
		cerr << "Output is: \"" << od.identifier << "\"" << endl;

		if (!plugin->initialise(channels, stepSize, blockSize)) {
			cerr << "ERROR: Plugin initialise (channels = " << channels
				<< ", stepSize = " << stepSize << ", blockSize = "
				<< blockSize << ") failed." << endl;
			goto done;
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for (int i = 0; i < parameters.size(); i++)			//MY OWN CODE 
		{
			if (parameters[i].identifier == "featureType");		// CHECKS IF THE PLUGIN HAS THE featureType permaeter and if it does sets it to the one i want
			plugin->setParameter("featureType", 3);
		}



		wrapper = dynamic_cast<PluginWrapper*>(plugin);
		if (wrapper) {
			// See documentation for
			// PluginInputDomainAdapter::getTimestampAdjustment
			PluginInputDomainAdapter* ida =
				wrapper->getWrapper<PluginInputDomainAdapter>();
			if (ida) adjustment = ida->getTimestampAdjustment();
		}

		// Here we iterate over the frames, avoiding asking the numframes in case it's streaming input.
		do {

			int count;

			if ((blockSize == stepSize) || (currentStep == 0)) {
				// read a full fresh block
				if ((count = sf_readf_float(sndfile, filebuf, blockSize)) < 0) {
					cerr << "ERROR: sf_readf_float failed: " << sf_strerror(sndfile) << endl;
					break;
				}
				if (count != blockSize) --finalStepsRemaining;
			}
			else {
				//  otherwise shunt the existing data down and read the remainder.
				memmove(filebuf, filebuf + (stepSize * channels), overlapSize * channels * sizeof(float));
				if ((count = sf_readf_float(sndfile, filebuf + (overlapSize * channels), stepSize)) < 0) {
					cerr << "ERROR: sf_readf_float failed: " << sf_strerror(sndfile) << endl;
					break;
				}
				if (count != stepSize) --finalStepsRemaining;
				count += overlapSize;
			}

			for (int c = 0; c < channels; ++c) {
				int j = 0;
				while (j < count) {
					plugbuf[c][j] = filebuf[j * sfinfo.channels + c];
					++j;
				}
				while (j < blockSize) {
					plugbuf[c][j] = 0.0f;
					++j;
				}
			}

			rt = RealTime::frame2RealTime(currentStep * stepSize, sfinfo.samplerate);
			
			features = plugin->process(plugbuf, rt);
			
			printFeatures
			(RealTime::realTime2Frame(rt + adjustment, sfinfo.samplerate),
				sfinfo.samplerate, od, outputNo, features, out, useFrames);
				
			if (sfinfo.frames > 0) {
				int pp = progress;
				progress = (int)((float(currentStep * stepSize) / sfinfo.frames) * 100.f + 0.5f);
				if (progress != pp && out) {
					cerr << "\r" << progress << "%";
				}
			}

			++currentStep;

		} while (finalStepsRemaining > 0);

		if (out) cerr << "\rDone" << endl;

		rt = RealTime::frame2RealTime(currentStep * stepSize, sfinfo.samplerate);

		features = plugin->getRemainingFeatures();
		
		printFeatures(RealTime::realTime2Frame(rt + adjustment, sfinfo.samplerate),
			sfinfo.samplerate, od, outputNo, features, out, useFrames);
			
		returnValue = 0;

	done:
		delete plugin;
		if (out) {
			out->close();
			delete out;
		}
		sf_close(sndfile);
		return returnValue;
	}


}

std::map <int, SegmentDetails> segmentDetails;
void fire()
{

	float totaltime = 0;
	int counter = 0;
	int segmentCount = 0;
	float lastTime =0;
	float lastSegmentTime =0;

	
	float temptime=0;
	float dt=0;
	Timer time;
	time.reset();
	
	bool doOnce = true;
		do {
			
			dt = (time.elapsed() - lastTime);
		
			if (combined.at(counter) > lastTime&& combined.at(counter) < ((totaltime)))
			{
				cout << "Fire ";
				counter++;
			}
			if (segmentCount >= combinedSegment.size())
			{
				if (doOnce == true)
				{
					cout << endl << segmentNames.back() << endl;
					doOnce = false;
				}
				
			}
			else if (combinedSegment.at(segmentCount) >= lastSegmentTime && combinedSegment.at(segmentCount) < ((totaltime)))
			{
				if (segmentDetails.at(segmentNames.at(segmentCount)).HighIntensity == true)
				{
					cout << endl << segmentNames.at(segmentCount) <<"   ------ PREDICTED CHORUS -------- "<< endl;
				}
				else
				{
					cout << endl << segmentNames.at(segmentCount) << endl;
				}
				segmentCount++;
			}

			lastSegmentTime = totaltime;
			lastTime = (totaltime);
			totaltime += dt;
			
		} while (totaltime < combined.back());
		

}




int main()
{ 
	 //cerr << PluginHostAdapter::getPluginPath().front();
	auto var = fs::current_path().string();
	system("set VAMP_PATH=. ;");

	cout << Vamp::PluginHostAdapter::getPluginPath()[0] << endl;

	string SongFileName;
	std::string path("");
	std::string ext(".wav");
	for (auto& p : fs::recursive_directory_iterator(path))
	{
		if (p.path().extension() == ext)
		{
			//SongFileName = p.path().string();
			SongFileName = fs::absolute(p).string();
			break;
		}
	}
	int seed;
	//GetFullPathName(SongFileName, MAX_PATH, )
	srand((unsigned)time(NULL));
	seed = rand() % 9;
	ofstream SettingsFile;
	SettingsFile.open("SettingsFile.txt");
	SettingsFile << SongFileName <<"\n";
	SettingsFile << seed;


	SettingsFile.close();



	run("test", "qm-vamp-plugins.dll", "qm-barbeattracker", "Beats", 0, SongFileName, "out.txt", false);
	Beats = false;
	Segments = true;
	run("test", "qm-vamp-plugins.dll", "qm-segmenter", "segmentation", 0, SongFileName, "Segments.txt", false);
	
	Segments = false;
	intensity = true;
	run("BBCTest", "bbc-vamp-plugins.dll", "bbc-intensity", "Intensity", 0, SongFileName, "intensity.txt", false);//Subconscious Grazer_Data/StreamingAssets/

	int size = miliSeconds.size();
	for (int i = 0; i < size; i++)
	{
		combined.push_back((seconds.back() + (miliSeconds.back() / 1000)));
		seconds.pop_back();
		miliSeconds.pop_back();
	}
	std::reverse(combined.begin(), combined.end());	

	size = miliSecondsSegment.size();

	for (int i = 0; i <size; i++)
	{
		combinedSegment.push_back((secondsSegment.back() + (miliSecondsSegment.back() / 1000)));
		secondsSegment.pop_back();
		miliSecondsSegment.pop_back();
	}
	std::reverse(combinedSegment.begin(), combinedSegment.end());

	size = miliSecondsIntensity.size();

	for (int i = 0; i <size; i++)
	{
		combinedIntensity.push_back((secondsIntensity.back() + (miliSecondsIntensity.back() / 1000)));
		secondsIntensity.pop_back();
		miliSecondsIntensity.pop_back();
	}
	std::reverse(combinedIntensity.begin(), combinedIntensity.end());


	for (int i = 0; i < combinedSegment.size(); i++)
	{
		cout << combinedSegment.at(i) << endl;
	}



	songDuration = combinedSegment.back() + lastDuration;

	/*
	for (all segments in segments)
		for (all intensity values in current segment)
			total += intensityValue
		next
		segmentAveragesVector.push(total)
		total = 0
	next
		*/
	float previousSegment = 0;
	float nextSegment = 0;
	float total = 0;
	int j = 0;

	std::vector<float> segmentIntensityMean;

	for (int i = 0; i < segmentNames.size(); i++)
	{
		previousSegment = combinedSegment.at(i);
		
		if ((i == (segmentNames.size()-1)))
		{
			//nextSegment = combinedSegment.back();
			nextSegment += 100000;
		}
		else
		{
			nextSegment = combinedSegment.at(i + 1);
		}
		int counter = 0;
		for (j; j < IntensityValues.size(); j++)
		{
			if (combinedIntensity.at(j) < nextSegment)
			{
				total += IntensityValues.at(j);
				counter++;
			}
			else
			{
				break;
			}
		}
		float mean = total / counter;
		total = 0;
		segmentIntensityMean.push_back(mean);
	}

	for (int i = 0; i < segmentIntensityMean.size(); i++)
	{
		cout << segmentNames.at(i) << " with intensity value: " << segmentIntensityMean.at(i) << endl;
	}

	
/*	SegmentDetails test;

	test.HighIntensity = true;
	test.Repeats = true;

	segmentDetails.insert(str_Seg(segmentNames.at(0), test));*/

	std::vector<float> segmentCalculations; // = segmentIntensityMean;
	int workingSegment = 0;
	int workingTotal = 0;
	int counter = 0;
	int max = 0;
	for (int j = 0; j < segmentNames.size(); j++)
	{
		if (max <= segmentNames.at(j))
		{
			max = segmentNames.at(j);
		}
	}
	for (int j = 1; j <= max; j++)
	{
		workingSegment = j;
		for (int i = 0; i < segmentIntensityMean.size(); i++)
		{
			if (segmentNames.at(i) == workingSegment)
			{
				workingTotal += segmentIntensityMean.at(i);
				counter++;
				
			}
		}
		SegmentDetails temp;
		temp.AverageIntensity = workingTotal / counter;
		if (counter > 1)
		{
			temp.repeats = true;
		}
		segmentDetails.insert(int_Seg(j, temp));
		counter = 0;
		workingTotal = 0;
	}

	std::vector<float> SortedValues;
	for (int i = 1; i <= segmentDetails.size(); i++)
	{
		SortedValues.push_back(segmentDetails.at(i).AverageIntensity);
	}

	std::sort(SortedValues.begin(), SortedValues.end());
	for (int j = 0; j < 2; j++)
	{
		for (int i = 1; i <= segmentDetails.size(); i++)
		{
			if (segmentDetails.at(i).AverageIntensity == SortedValues.back())
			{
				segmentDetails.at(i).HighIntensity = true;
			
			}
		}
		SortedValues.pop_back();
	}

	fire();
	
	
}

///////
/*

Chorus - Repeated, High Energy

Verse - Repeated Low Energy.


First, check sections for Repeats.

	High Value Repeats Are Choruses/Drops/Post Chorus. 

	Low are verses

Second, 
	If two High Energy Are next to each other then one is post chorus.

Third,
	If one before Chorus is lower than Chorus but higher than verse then its a build

Extras are special sections/Verses/Intro/Outro

/////////////////////////////////////////////////////


Average Values between Same Segments.
Sort Vector From highest to lowest.
	This is now a vector that shows the intensity of each type of segment.

Highest values are Choruses



*/