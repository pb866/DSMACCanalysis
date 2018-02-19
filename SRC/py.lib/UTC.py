from timezonefinder import TimezoneFinder
from pytz import timezone
import pytz
from datetime import datetime

tf = TimezoneFinder()

def offset(latitude, longitude):
    today = datetime.now()
    tz_target = timezone(tf.certain_timezone_at(lat=latitude, lng=longitude))
    # ATTENTION: tz_target could be None! handle error case
    today_target = tz_target.localize(today)
    today_utc = pytz.utc.localize(today)
    return (today_utc - today_target).total_seconds()
