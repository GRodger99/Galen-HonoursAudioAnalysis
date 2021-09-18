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
struct SegmentDetails	//a nice way to store the details about each segment
{
	bool repeats = false;
	bool HighIntensity = false;
	float AverageIntensity = 0;
};
typedef std::pair<int, SegmentDetails> int_Seg;	//give each segment type an identifier


#define HOST_VERSION "1"


std::vector<float> combined;	//create all the containers needed
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

			if (out && Beats == true)	//swap depending on what plugins are being run
			{
				seconds.push_back(rt.sec);	//the values are stored seperately so i am needed to store them seperately and combine them later
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

			if (f.hasDuration) {
				rt = f.duration;
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

		//change the vamp parameters of the plugins being used.
		for (int i = 0; i < parameters.size(); i++)	
		{
			if (parameters[i].identifier == "featureType");		// checks if the plugin has featureType permaeter and if it does sets it to 3 which is to use Tibral Segments
			plugin->setParameter("featureType", 3);

			if (parameters[i].identifier == "neighbourhoodLimit");		//then check for the next parameter needed. this one controls minimum segment length
			plugin->setParameter("neighbourhoodLimit", 12);
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

int main()
{ 
	
	_putenv("VAMP_PATH=./Plugins");	//create an environment variable for vamp to know where the plugins are stored
	string SongFileName;
	std::string path("");
	std::string ext(".wav");
	for (auto& p : fs::recursive_directory_iterator(path))	//search through evertyhing in the current folder andcheck if it is a wav file
	{
		if (p.path().extension() == ext)
		{
			SongFileName = fs::absolute(p).string();	//if it is the return the name of the file 
			break;//stop after finding the first wav
		}
	}



	run("test", "qm-vamp-plugins.dll", "qm-barbeattracker", "Beats", 0, SongFileName, "Subconscious Grazer_Data/StreamingAssets/out.txt", false);	//run the QM beat tracking plugin 
	Beats = false;	//set the beats to false and set segments to true so that the segments are poroperly stored
	Segments = true;
	run("test", "qm-vamp-plugins.dll", "qm-segmenter", "segmentation", 0, SongFileName, "Subconscious Grazer_Data/StreamingAssets/Segments.txt", false);	//run segments
	
	Segments = false;	//setup for intensity
	intensity = true;
	run("BBCTest", "bbc-vamp-plugins.dll", "bbc-intensity", "Intensity", 0, SongFileName, "Subconscious Grazer_Data/StreamingAssets/intensity.txt", false);	//then finally run the intensity 

	int size = miliSeconds.size();	//get the total number of beats
	for (int i = 0; i < size; i++) //for each beat
	{
		combined.push_back((seconds.back() + (miliSeconds.back() / 1000)));	//store the values from the individual vectors and merge them into one
		seconds.pop_back();
		miliSeconds.pop_back();
	}
	std::reverse(combined.begin(), combined.end());		//the re-order it so that the first beat is values 0

	size = miliSecondsSegment.size(); //do the same thing but for the segments instead

	for (int i = 0; i <size; i++)
	{
		combinedSegment.push_back((secondsSegment.back() + (miliSecondsSegment.back() / 1000)));	//merge the segments
		secondsSegment.pop_back();
		miliSecondsSegment.pop_back();
	}
	std::reverse(combinedSegment.begin(), combinedSegment.end());

	size = miliSecondsIntensity.size();	//then for the intensity

	for (int i = 0; i <size; i++)
	{
		combinedIntensity.push_back((secondsIntensity.back() + (miliSecondsIntensity.back() / 1000)));
		secondsIntensity.pop_back();
		miliSecondsIntensity.pop_back();
	}
	std::reverse(combinedIntensity.begin(), combinedIntensity.end());	//keep it consistant with the others

	float previousSegment = 0;
	float nextSegment = 0;
	float total = 0;
	int j = 0;

	std::vector<float> segmentIntensityMean;	//create the storage for the average intensity value of all the segments

	for (int i = 0; i < segmentNames.size(); i++)	//loop through every segment
	{
		previousSegment = combinedSegment.at(i);
		
		if ((i == (segmentNames.size()-1)))// if it is the last segment
		{
			nextSegment += 100000;	//make it some large value to make sure all the intensity values are used
		}
		else
		{
			nextSegment = combinedSegment.at(i + 1);	//set the next segment to be i+1
		}
		int counter = 0;
		for (j; j < IntensityValues.size(); j++)	//after setting up the next segments timer
		{
			if (combinedIntensity.at(j) < nextSegment)	//if the current intensity time is before the next segment
			{
				total += IntensityValues.at(j);	//then add the value to that segment's total
				counter++;
			}
			else
			{
				break;
			}
		}
		float mean = total / counter;	//calculate the mean
		total = 0;	//the reset the total
		segmentIntensityMean.push_back(mean);	//add this new mean to the vector of mean values that is in the same order as the segments
	}

	std::vector<float> segmentCalculations; 
	int workingSegment = 0;
	int workingTotal = 0;
	int counter = 0;
	int max = 0;
	for (int j = 0; j < segmentNames.size(); j++)//go through every segment and find out the total number of different segments there are
	{
		if (max <= segmentNames.at(j))
		{
			max = segmentNames.at(j);
		}
	}
	for (int j = 1; j <= max; j++)	//go through every type of segment
	{
		workingSegment = j;
		for (int i = 0; i < segmentIntensityMean.size(); i++)	//get the average value for every type of this segment.
		{
			if (segmentNames.at(i) == workingSegment)	//loop through all the segments and see if there are multiple of this one
			{
				workingTotal += segmentIntensityMean.at(i);	// if there are then add them
				counter++;
				
			}
		}
		SegmentDetails temp;
		temp.AverageIntensity = workingTotal / counter;	//get the average intensity for this type of segment

		if (counter > 1)	//if it calculated the mean for this type of segment more than once
		{
			temp.repeats = true;	
		}
		segmentDetails.insert(int_Seg(j, temp));	//add it to the map using the segment number as the index
		counter = 0;
		workingTotal = 0;
	}
	//Calulating the Chorus using intensity
	std::vector<float> SortedValues;
	for (int i = 1; i <= segmentDetails.size(); i++)
	{
		SortedValues.push_back(segmentDetails.at(i).AverageIntensity);
	}

	std::sort(SortedValues.begin(), SortedValues.end());	
	for (int j = 0; j < 2; j++)	//loop through and take only 2 values
	{
		for (int i = 1; i <= segmentDetails.size(); i++)
		{
			if (segmentDetails.at(i).AverageIntensity == SortedValues.back())	//if the average intensity of the current segment is the same as the highest value
			{
				segmentDetails.at(i).HighIntensity = true; //then set it to be a predicted chorus
			
			}
		}
		SortedValues.pop_back();	//remove the value from the back
	}

	//outputting for the "settingsFile" which has the specific values
	vector<float> Outputted;
	int outputCount=0;
	int seed;	
	srand((unsigned)time(NULL));
	seed = rand() % 9;	//create the seed in this program so that it is the same whenever you use the same data
	ofstream SettingsFile;
	SettingsFile.open("Subconscious Grazer_Data/StreamingAssets/SettingsFile.txt");
	SettingsFile << SongFileName << "\n";	//add the song path to the file
	SettingsFile << seed <<endl;	//also add the seed

	for (int i = 0; i < segmentNames.size(); i++)	
	{
		if (segmentDetails.at(segmentNames.at(i)).HighIntensity == true)	//if the first type of segment is a predicted chorus
		{
			if (outputCount < 1 )	//only do this for the first type to output
			{
				SettingsFile << segmentNames.at(i) << endl;	//output the name if that segment
				Outputted.push_back(segmentNames.at(i));	//add it to a vector of the type fo segments that have been outputted
				outputCount++;
			}
			else
			{
				if (Outputted.at(0) == segmentNames.at(i))	//if the first value outputted is the same as the current one
				{
					break;	//dont do anything
				}
				else
				{
					SettingsFile << segmentNames.at(i) << endl;	//output the second one 
					Outputted.push_back(segmentNames.at(i));
					outputCount++;
				}
			}
		}
	}

	SettingsFile.close();	
}
