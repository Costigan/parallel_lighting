namespace viper.corelib.horizon
{
    public enum ProcessingQueueType { Horizon = 0, AverageSun = 1, AverageEarth = 2, AverageSunEarth = 3 };

    public class ProcessingQueue
    {
        public string Name { get; set; }
        public ProcessingQueueType Type { get; set; } = ProcessingQueueType.Horizon;
      
        // Horizon Queue
    }
}
